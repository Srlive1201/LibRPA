#include "parallel_mpi.h"
#include "profiler.h"
#include "matrix_m_parallel_utils.h"
#include "params.h"
#include "constants.h"
#include "exx.h"
#include "pbc.h"
#include "lapack_connector.h"
#include "vector3_order.h"
#include "libri_utils.h"
#include "stl_io_helper.h"
#ifdef LIBRPA_USE_LIBRI
#include <RI/physics/Exx.h>
// #include "print_stl.h"
#endif

namespace LIBRPA
{

ComplexMatrix Exx::get_dmat_cplx_R_global(const int& ispin, const Vector3_Order<int>& R)
{
    const auto nspins = this->mf_.get_n_spins();
    auto dmat_cplx = this->mf_.get_dmat_cplx_R(ispin, this->kfrac_list_, R);
    // renormalize to single spin channel
    dmat_cplx *= 0.5 * nspins;
    return dmat_cplx;
}

ComplexMatrix Exx::extract_dmat_cplx_R_IJblock(const ComplexMatrix& dmat_cplx, const atom_t& I, const atom_t& J)
{
    const auto I_num = atom_nw.at(I);
    const auto J_num = atom_nw.at(J);
    ComplexMatrix dmat_cplx_IJR(I_num, J_num);
    for (size_t i = 0; i != I_num; i++)
    {
        size_t i_glo = atom_iw_loc2glo(I, i);
        for (size_t j = 0; j != J_num; j++)
        {
            size_t j_glo = atom_iw_loc2glo(J, j);
            dmat_cplx_IJR(i, j) = dmat_cplx(i_glo, j_glo);
        }
    }
    return dmat_cplx_IJR;
}

void Exx::build_dmat_R(const Vector3_Order<int> &R)
{
    const auto nspins = this->mf_.get_n_spins();

    for (int is = 0; is != nspins; is++)
    {
        auto dmat_cplx = this->get_dmat_cplx_R_global(is, R);
        for ( int I = 0; I != natom; I++)
        {
            for (int J = 0; J != natom; J++)
            {
                const auto dmat_cplx_IJR = this->extract_dmat_cplx_R_IJblock(dmat_cplx, I, J);
                this->warn_dmat_IJR_nonzero_imag(dmat_cplx_IJR, is, I, J, R);
                this->dmat[is][I][J][R] = std::make_shared<matrix>();
                *(this->dmat[is][I][J][R]) = dmat_cplx_IJR.real();
            }
        }
    }
}

void Exx::build_dmat_R(const atom_t& I, const atom_t& J, const Vector3_Order<int> &R)
{
    const auto nspins = this->mf_.get_n_spins();

    for (int is = 0; is != nspins; is++)
    {
        auto dmat_cplx = this->get_dmat_cplx_R_global(is, R);
        const auto dmat_cplx_IJR = this->extract_dmat_cplx_R_IJblock(dmat_cplx, I, J);
        this->warn_dmat_IJR_nonzero_imag(dmat_cplx_IJR, is, I, J, R);
        this->dmat[is][I][J][R] = std::make_shared<matrix>();
        *(this->dmat[is][I][J][R]) = dmat_cplx_IJR.real();
    }
}

void Exx::warn_dmat_IJR_nonzero_imag(const ComplexMatrix& dmat_cplx, const int& ispin, const atom_t& I, const atom_t& J, const Vector3_Order<int> R)
{
    if (dmat_cplx.get_max_abs_imag() > 1e-2)
        printf("Warning: complex-valued density matrix, spin %d IJR %zu %zu (%d, %d, %d)\n", ispin, I, J, R.x, R.y, R.z);
}

void Exx::build_exx_orbital_energy(const atpair_R_mat_t &LRI_Cs,
                                   const vector<Vector3_Order<int>> &Rlist,
                                   const Vector3_Order<int> &R_period,
                                   const atpair_R_mat_t &coul_mat)
{
    if (Params::use_libri_exx)
        this->build_exx_orbital_energy_LibRI(LRI_Cs, Rlist, R_period, coul_mat);
    else
    {
        if (mpi_comm_world_h.is_root())
            printf("Error: EXX orbital energy without LibRI is not implemented yet\n");
        mpi_comm_world_h.barrier();
        throw std::logic_error("not implemented");
    }
}

void Exx::build_exx_orbital_energy_LibRI(const atpair_R_mat_t &LRI_Cs,
                                         const vector<Vector3_Order<int>> &Rlist,
                                         const Vector3_Order<int> &R_period,
                                         const atpair_R_mat_t &coul_mat)
{
    const auto& n_aos = this->mf_.get_n_aos();
    const auto& n_spins = this->mf_.get_n_spins();
    const auto& n_kpts = this->mf_.get_n_kpoints();
    const auto& n_bands = this->mf_.get_n_bands();

#ifdef LIBRPA_USE_LIBRI
    if (mpi_comm_world_h.is_root())
        printf("Computing EXX orbital energy using LibRI\n");
    mpi_comm_world_h.barrier();
    RI::Exx<int, int, 3, double> exx_libri;
    map<int,std::array<double,3>> atoms_pos;
    for(int i=0;i!=atom_mu.size();i++)
        atoms_pos.insert(pair<int,std::array<double,3>>{i,{0,0,0}});
    std::array<double,3> xa{latvec.e11,latvec.e12,latvec.e13};
    std::array<double,3> ya{latvec.e21,latvec.e22,latvec.e23};
    std::array<double,3> za{latvec.e31,latvec.e32,latvec.e33};
    std::array<std::array<double,3>,3> lat_array{xa,ya,za};
    std::array<int,3> period_array{R_period.x,R_period.y,R_period.z};
    exx_libri.set_parallel(MPI_COMM_WORLD, atoms_pos, lat_array, period_array);
    exx_libri.set_csm_threshold(Params::libri_exx_threshold_CSM);

    // Initialize Cs libRI container on each process
    // Note: we use different treatment in different routings
    //     R-tau routing:
    //         Each process has a full Cs copy.
    //         Thus in each process we only pass a few to LibRI container.
    //     atom-pair routing:
    //         Cs is already distributed across all processes.
    //         Pass the all Cs to libRI container.
    std::map<int, std::map<std::pair<int,std::array<int,3>>,RI::Tensor<double>>> Cs_libri;

    Profiler::start("build_exx_orbital_energy_1", "Prepare C libRI object");
    if(exx_parallel_type == parallel_type::LIBRI_USED)
    {
        if (chi_parallel_type == parallel_type::R_TAU)
        {
            std::vector<std::pair<atom_t, std::pair<atom_t, Vector3_Order<int>>>> Cs_IJRs_local;
            size_t n_Cs_IJRs = 0;
            size_t n_Cs_IJRs_local = 0;

            for(auto &Ip:LRI_Cs)
                for(auto &Jp:Ip.second)
                    for(auto &Rp:Jp.second)
                    {
                        const auto &I = Ip.first;
                        const auto &J = Jp.first;
                        const auto &R=Rp.first;
                        if ((n_Cs_IJRs++) % mpi_comm_world_h.nprocs == mpi_comm_world_h.myid)
                        {
                            Cs_IJRs_local.push_back({I, {J, R}});
                            n_Cs_IJRs_local++;
                        }
                    }
            printf("| Number of Cs on Proc %4d: %zu\n", mpi_comm_world_h.myid, n_Cs_IJRs_local);
            for (auto &IJR: Cs_IJRs_local)
            {
                const auto I = IJR.first;
                const auto J = IJR.second.first;
                const auto R = IJR.second.second;
                const std::array<int,3> Ra{R.x,R.y,R.z};
                const auto mat = transpose(*(LRI_Cs.at(I).at(J).at(R)));
                std::valarray<double> mat_array(mat.c, mat.size);
                std::shared_ptr<std::valarray<double>> mat_ptr = std::make_shared<std::valarray<double>>();
                *mat_ptr=mat_array;
                Cs_libri[I][{J, Ra}] = RI::Tensor<double>({atom_mu[I], atom_nw[I], atom_nw[J]}, mat_ptr);
            }
        }
        else if (chi_parallel_type == parallel_type::ATOM_PAIR)
        {
            for (auto &I_JRCs: LRI_Cs)
            {
                const auto &I = I_JRCs.first;
                for (auto &J_RCs: I_JRCs.second)
                {
                    const auto &J = J_RCs.first;
                    for (auto &R_Cs: J_RCs.second)
                    {
                        const auto &R = R_Cs.first;
                        const std::array<int,3> Ra{R.x,R.y,R.z};
                        const auto mat = transpose(*R_Cs.second);
                        std::valarray<double> mat_array(mat.c, mat.size);
                        std::shared_ptr<std::valarray<double>> mat_ptr = std::make_shared<std::valarray<double>>();
                        *mat_ptr=mat_array;
                        Cs_libri[I][{J, Ra}] = RI::Tensor<double>({atom_mu[I], atom_nw[I], atom_nw[J]}, mat_ptr);
                    }
                }
            }
        }
    }
    else
    {
        // Not implemented
        MPI_Wrapper::barrier_world();
        throw std::logic_error("not implemented");
    }
    fout_para << "Number of Cs keys: " << get_num_keys(Cs_libri) << "\n";
    print_keys(LIBRPA::fout_para, Cs_libri);
    exx_libri.set_Cs(Cs_libri, Params::libri_exx_threshold_C);
    Profiler::stop("build_exx_orbital_energy_1");
    printf("Task %4d: C setup for EXX\n", mpi_comm_world_h.myid);

    // initialize Coulomb matrix
    Profiler::start("build_exx_orbital_energy_2", "Prepare V libRI object");
    std::map<int, std::map<std::pair<int,std::array<int,3>>,RI::Tensor<double>>> V_libri;
    for (auto IJR: dispatch_vector_prod(get_atom_pair(coul_mat), Rlist, mpi_comm_world_h.myid, mpi_comm_world_h.nprocs, true, true))
    {
        const auto I = IJR.first.first;
        const auto J = IJR.first.second;
        const auto R = IJR.second;
        const auto& VIJR = coul_mat.at(I).at(J).at(R);
        std::array<int,3> Ra{R.x,R.y,R.z};
        std::valarray<double> VIJR_va(VIJR->c, VIJR->size);
        auto pv = std::make_shared<std::valarray<double>>();
        *pv = VIJR_va;
        V_libri[I][{J, Ra}] = RI::Tensor<double>({size_t(VIJR->nr), size_t(VIJR->nc)}, pv);
    }
    fout_para << "Number of V keys: " << get_num_keys(V_libri) << "\n";
    exx_libri.set_Vs(V_libri, Params::libri_exx_threshold_V);
    Profiler::stop("build_exx_orbital_energy_2");
    printf("Task %4d: V setup for EXX\n", mpi_comm_world_h.myid);
    // cout << V_libri << endl;

    // initialize density matrix
    vector<atpair_t> atpair_dmat;
    for (int I = 0; I < atom_nw.size(); I++)
        for (int J = 0; J < atom_nw.size(); J++)
            atpair_dmat.push_back({I, J});
    const auto dmat_IJRs_local = dispatch_vector_prod(atpair_dmat, Rlist, mpi_comm_world_h.myid, mpi_comm_world_h.nprocs, true, true);

    // prepare scalapack array descriptors
    Array_Desc desc_nao_nao(LIBRPA::blacs_ctxt_world_h);
    Array_Desc desc_nband_nao(LIBRPA::blacs_ctxt_world_h);
    Array_Desc desc_nband_nband(LIBRPA::blacs_ctxt_world_h);
    Array_Desc desc_nband_nband_fb(LIBRPA::blacs_ctxt_world_h);
    desc_nao_nao.init_1b1p(n_aos, n_aos, 0, 0);
    desc_nband_nao.init_1b1p(n_bands, n_aos, 0, 0);
    desc_nband_nband.init_1b1p(n_bands, n_bands, 0, 0);
    desc_nband_nband_fb.init(n_bands, n_bands, n_bands, n_bands, 0, 0);
    // local 2D-block submatrices
    auto Hexx_nao_nao = init_local_mat<complex<double>>(desc_nao_nao, MAJOR::COL);
    auto temp_nband_nao = init_local_mat<complex<double>>(desc_nband_nao, MAJOR::COL);
    auto Hexx_nband_nband = init_local_mat<complex<double>>(desc_nband_nband, MAJOR::COL);
    auto Hexx_nband_nband_fb = init_local_mat<complex<double>>(desc_nband_nband_fb, MAJOR::COL);
    // printf("%d size naonao %d nbandnao %d nbandnband %d\n",
    //        blacs_ctxt_world_h.myid, Hexx_nao_nao.size(), temp_nband_nao.size(), Hexx_nband_nband.size());

    for (auto isp = 0; isp != n_spins; isp++)
    {
        Profiler::start("build_exx_orbital_energy_3", "Prepare DM libRI object");
        std::map<int, std::map<std::pair<int,std::array<int,3>>,RI::Tensor<double>>> dmat_libri;
        for (const auto &R: Rlist)
        {
            std::array<int,3> Ra{R.x,R.y,R.z};
            const auto dmat_cplx = this->get_dmat_cplx_R_global(isp, R);
            for (const auto& IJR: dmat_IJRs_local)
            {
                if (IJR.second == R)
                {
                    const auto& I = IJR.first.first;
                    const auto& J = IJR.first.second;
                    const auto dmat_IJR = this->extract_dmat_cplx_R_IJblock(dmat_cplx, I, J);
                    this->warn_dmat_IJR_nonzero_imag(dmat_IJR, isp, I, J, R);
                    std::valarray<double> dmat_va(dmat_IJR.real().c, dmat_IJR.size);
                    auto pdmat = std::make_shared<std::valarray<double>>();
                    *pdmat = dmat_va;
                    dmat_libri[I][{J, Ra}] = RI::Tensor<double>({size_t(dmat_IJR.nr), size_t(dmat_IJR.nc)}, pdmat);
                }
            }
        }
        fout_para << "Number of Dmat keys: " << get_num_keys(dmat_libri) << "\n";
        print_keys(LIBRPA::fout_para, dmat_libri);
        exx_libri.set_Ds(dmat_libri, Params::libri_exx_threshold_D);
        Profiler::stop("build_exx_orbital_energy_3");
        printf("Task %4d: DM setup for EXX\n", mpi_comm_world_h.myid);

        Profiler::start("build_exx_orbital_energy_4", "Call libRI Hexx calculation");
        exx_libri.cal_Hs();
        Profiler::stop("build_exx_orbital_energy_4");
        printf("Task %4d: cal_Hs elapsed time: %f\n", mpi_comm_world_h.myid, Profiler::get_wall_time_last("build_exx_orbital_energy_4"));
        print_keys(LIBRPA::fout_para, exx_libri.Hs);
        LIBRPA::fout_para << "exx_libri.Hs:\n" << exx_libri.Hs << endl;

        // collect necessary data
        Profiler::start("build_exx_orbital_energy_5", "Collect Hexx IJ from world");
        // collect the IJ pair of Hs with all R, do the Fourier transform
        const auto set_IJ_naonao = get_necessary_IJ_from_block_2D(atomic_basis_wfc,
                                                                  atomic_basis_wfc,
                                                                  desc_nao_nao);
        const auto Iset_Jset = convert_IJset_to_Iset_Jset(set_IJ_naonao);
        const auto I_JallR_Hs = RI::Communicate_Tensors_Map_Judge::comm_map2_first(LIBRPA::mpi_comm_world_h.comm,
                exx_libri.Hs, Iset_Jset.first, Iset_Jset.second);
        Profiler::stop("build_exx_orbital_energy_5");
        printf("Task %4d: tensor communicate elapsed time: %f\n", mpi_comm_world_h.myid, Profiler::get_wall_time_last("build_exx_orbital_energy_5"));
        // cout << I_JallR_Hs << endl;

        for (int ik = 0; ik < n_kpts; ik++)
        {
            Hexx_nao_nao.zero_out();
            Profiler::start("build_exx_orbital_energy_6", "Hexx IJ -> 2D block");
            const Vector3_Order<double>& kfrac = this->kfrac_list_[ik];
            for (auto &T: I_JallR_Hs)
            {
                for (auto &T1: T.second)
                {
                    const auto& I = T.first;
                    const auto& J = T1.first.first;
                    const auto& Ra = T1.first.second;
                    double ang = kfrac * Vector3_Order<int>{Ra[0], Ra[1], Ra[2]} * TWO_PI;
                    complex<double> kphase = complex<double>(cos(ang), sin(ang));
                    const auto alpha = kphase / double(n_kpts);
                    collect_block_from_IJ_storage(Hexx_nao_nao, desc_nao_nao,
                                                  atomic_basis_wfc, atomic_basis_wfc,
                                                  I, J, alpha, T1.second.ptr(), MAJOR::ROW);
                }
            }
            Profiler::stop("build_exx_orbital_energy_6");
            // printf("%s\n", str(Hexx_nao_nao).c_str());
            const auto &wfc_isp_k = this->mf_.get_eigenvectors()[isp][ik];
            blacs_ctxt_world_h.barrier();
            const auto wfc_block = get_local_mat(wfc_isp_k.c, MAJOR::ROW, desc_nband_nao, MAJOR::COL).conj();
            // printf("%s\n", str(wfc_block).c_str());
            // printf("%s\n", desc_nao_nao.info_desc().c_str());
            // printf("%s\n", desc_nband_nao.info_desc().c_str());
            Profiler::start("build_exx_orbital_energy_7", "Rotate Hexx ij -> KS");
            ScalapackConnector::pgemm_f('N', 'N', n_bands, n_aos, n_aos, 1.0,
                                        wfc_block.ptr(), 1, 1, desc_nband_nao.desc,
                                        Hexx_nao_nao.ptr(), 1, 1, desc_nao_nao.desc,
                                        0.0,
                                        temp_nband_nao.ptr(), 1, 1, desc_nband_nao.desc);
            ScalapackConnector::pgemm_f('N', 'C', n_bands, n_bands, n_aos, -1.0,
                                        temp_nband_nao.ptr(), 1, 1, desc_nao_nao.desc,
                                        wfc_block.ptr(), 1, 1, desc_nband_nao.desc,
                                        0.0,
                                        Hexx_nband_nband.ptr(), 1, 1, desc_nband_nband.desc);
            Profiler::stop("build_exx_orbital_energy_7");

            // collect to master
            Profiler::start("build_exx_orbital_energy_8", "Collect Eexx to root process");
            ScalapackConnector::pgemr2d_f(n_bands, n_bands,
                                          Hexx_nband_nband.ptr(), 1, 1, desc_nband_nband.desc,
                                          Hexx_nband_nband_fb.ptr(), 1, 1, desc_nband_nband_fb.desc,
                                          desc_nband_nband_fb.ictxt());
            // cout << "Hexx_nband_nband_fb isp " << isp  << " ik " << ik << endl << Hexx_nband_nband_fb;
            if (LIBRPA::blacs_ctxt_world_h.myid == 0)
                for (int ib = 0; ib != n_bands; ib++)
                    this->Eexx[isp][ik][ib] = Hexx_nband_nband_fb(ib, ib).real();
            Profiler::stop("build_exx_orbital_energy_8");
        }
    }
    // debug, print the Hexx matrices
    // for (const auto& isp_IJkH: this->Hexx)
    // {
    //     const auto& isp = isp_IJkH.first;
    //     for (const auto& I_JkH: isp_IJkH.second)
    //     {
    //         const auto& I = I_JkH.first;
    //         for (const auto& J_kH: I_JkH.second)
    //         {
    //             const auto& J = J_kH.first;
    //             for (const auto& k_H: J_kH.second)
    //             {
    //                 cout << isp << " " << I << " " << J << " {" << k_H.first << "} whole size: " << k_H.second->size << endl;
    //                 print_complex_matrix("", *k_H.second);
    //             }
    //         }
    //     }
    // }

    // Rotate the Hamiltonian in basis of KS state for each spin and kpoints
    // for (int isp = 0; isp < n_spins; isp++)
    // {
    //     for (int ik = 0; ik < n_kpts; ik++)
    //     {
    //         const auto& k = this->kfrac_list_[ik];
    //         const auto& IJHs = this->Hexx.at(isp).at(k);
    //         // retrieve the global Hexx
    //         ComplexMatrix hexx(n_aos, n_aos);
    //         for (const auto& I_JH: IJHs)
    //         {
    //             const auto& I = I_JH.first;
    //             const auto I_num = atom_nw[I];
    //             for (const auto& J_H: I_JH.second)
    //             {
    //                 const auto& J = J_H.first;
    //                 const auto J_num = atom_nw[J];
    //                 const auto& H = J_H.second;
    //                 for (int i = 0; i != I_num; i++)
    //                 {
    //                     size_t i_glo = atom_iw_loc2glo(I, i);
    //                     for (int j = 0; j != J_num; j++)
    //                     {
    //                         size_t j_glo = atom_iw_loc2glo(J, j);
    //                         hexx(i_glo, j_glo) = (*H)(i, j);
    //                     }
    //                 }
    //             }
    //         }
    //         const auto& wfc = this->mf_.get_eigenvectors()[isp][ik];
    //         this->Hexx_KS[isp][ik] = wfc * hexx * transpose(wfc, true);
    //         this->Hexx_KS[isp][ik] *= -1.0;
    //         for (int ib = 0; ib != this->mf_.get_n_bands(); ib++)
    //             this->Eexx[isp][ik][ib] = this->Hexx_KS[isp][ik](ib, ib).real();
    //     }
    // }
#else
    if (mpi_comm_world_h.is_root())
    {
        printf("Error: trying build EXX orbital energy with LibRI, but the program is not compiled against LibRI\n");
    }
    throw std::logic_error("compilation");
    mpi_comm_world_h.barrier();
#endif
}
}

