#include "exx.h"

#include "envs_io.h"
#include "envs_mpi.h"
#include "envs_blacs.h"
#include "utils_blacs.h"
#include "profiler.h"
#include "matrix_m_parallel_utils.h"
#include "params.h"
#include "constants.h"
#include "pbc.h"
#include "geometry.h"
#include "lapack_connector.h"
#include "vector3_order.h"
#include "libri_utils.h"
#include "stl_io_helper.h"
#ifdef LIBRPA_USE_LIBRI
#include <RI/physics/Exx.h>
#include <RI/ri/Cell_Nearest.h>
#else
#include "libri_stub.h"
#endif
#include "utils_io.h"

namespace LIBRPA
{

Exx::Exx(const MeanField &mf,
         const vector<Vector3_Order<double>> &kfrac_list,
         const Vector3_Order<int> &period)
    : mf_(mf), kfrac_list_(kfrac_list), period_(period)
{
    is_rspace_build_ = false;
    is_kspace_built_ = false;
};

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
        utils::lib_printf("Warning: complex-valued density matrix, spin %d IJR %zu %zu (%d, %d, %d)\n", ispin, I, J, R.x, R.y, R.z);
}


void Exx::build(const Cs_LRI &Cs,
                const vector<Vector3_Order<int>> &Rlist,
                const atpair_R_mat_t &coul_mat)
{
    using LIBRPA::envs::mpi_comm_global;
    using LIBRPA::envs::mpi_comm_global_h;

    assert (parallel_routing == ParallelRouting::LIBRI);

    if (this->is_rspace_build_)
    {
        return;
    }

    const auto& n_spins = this->mf_.get_n_spins();

#ifdef LIBRPA_USE_LIBRI
    if (mpi_comm_global_h.is_root())
    {
        utils::lib_printf("Computing EXX orbital energy using LibRI\n");
    }
    mpi_comm_global_h.barrier();

    RI::Exx<int, int, 3, double> exx_libri;
    map<int,std::array<double,3>> atoms_pos;
    for(int i=0;i!=atom_mu.size();i++)
        atoms_pos.insert(pair<int,std::array<double,3>>{i,{0,0,0}});

    std::array<double,3> xa{latvec.e11,latvec.e12,latvec.e13};
    std::array<double,3> ya{latvec.e21,latvec.e22,latvec.e23};
    std::array<double,3> za{latvec.e31,latvec.e32,latvec.e33};
    std::array<std::array<double,3>,3> lat_array{xa,ya,za};
    std::array<int,3> period_array{period_.x,period_.y,period_.z};
    exx_libri.set_parallel(mpi_comm_global, atoms_pos, lat_array, period_array);

    // Initialize Cs libRI container on each process
    // Note: we use different treatment in different routings
    //     R-tau routing:
    //         Each process has a full Cs copy.
    //         Thus in each process we only pass a few to LibRI container.
    //     atom-pair routing:
    //         Cs is already distributed across all processes.
    //         Pass the all Cs to libRI container.

    Profiler::start("build_real_space_exx_1", "Prepare C libRI object");
    envs::ofs_myid << "Number of Cs keys: " << get_num_keys(Cs.data_libri) << "\n";
    // print_keys(envs::ofs_myid, Cs.data_libri);
    exx_libri.set_Cs(Cs.data_libri, Params::libri_exx_threshold_C);
    Profiler::stop("build_real_space_exx_1");
    envs::ofs_myid << "Finished setup Cs for EXX\n";
    std::flush(envs::ofs_myid);

    // initialize Coulomb matrix
    Profiler::start("build_real_space_exx_2", "Prepare V libRI object");
    std::map<int, std::map<std::pair<int,std::array<int,3>>, RI::Tensor<double>>> V_libri;
    if (LIBRPA::parallel_routing == LIBRPA::ParallelRouting::R_TAU)
    {
        // Full Coulomb case, have to re-distribute
        for (auto IJR: dispatch_vector_prod(get_atom_pair(coul_mat), Rlist, mpi_comm_global_h.myid, mpi_comm_global_h.nprocs, true, true))
        {
            const auto I = IJR.first.first;
            const auto J = IJR.first.second;
            const auto R = IJR.second;
            const auto& VIJR = coul_mat.at(I).at(J).at(R);
            // debug
            // printf("I J R %zu %zu %d %d %d, max(V) %f\n", I, J, R.x, R.y, R.z, VIJR->max());
            std::array<int,3> Ra{R.x,R.y,R.z};
            std::valarray<double> VIJR_va(VIJR->c, VIJR->size);
            auto pv = std::make_shared<std::valarray<double>>();
            *pv = VIJR_va;
            V_libri[I][{J, Ra}] = RI::Tensor<double>({size_t(VIJR->nr), size_t(VIJR->nc)}, pv);
        }
    }
    else
    {
        for (const auto &I_JRV: coul_mat)
        {
            const auto I = I_JRV.first;
            for (const auto &J_RV: I_JRV.second)
            {
                const auto J = J_RV.first;
                for (const auto &R_V: J_RV.second)
                {
                    const auto &R = R_V.first;
                    const auto &V = R_V.second;
                    std::array<int,3> Ra{R.x,R.y,R.z};
                    std::valarray<double> VIJR_va(V->c, V->size);
                    auto pv = std::make_shared<std::valarray<double>>();
                    *pv = VIJR_va;
                    V_libri[I][{J, Ra}] = RI::Tensor<double>({size_t(V->nr), size_t(V->nc)}, pv);
                }
            }
        }
    }
    envs::ofs_myid << "Number of V keys: " << get_num_keys(V_libri) << "\n";
    exx_libri.set_Vs(V_libri, Params::libri_exx_threshold_V);
    Profiler::stop("build_real_space_exx_2");
    utils::lib_printf("Task %4d: V setup for EXX\n", mpi_comm_global_h.myid);
    // cout << V_libri << endl;

    // initialize density matrix
    vector<atpair_t> atpair_dmat;
    for (int I = 0; I < atom_nw.size(); I++)
        for (int J = 0; J < atom_nw.size(); J++)
            atpair_dmat.push_back({I, J});
    const auto dmat_IJRs_local = dispatch_vector_prod(atpair_dmat, Rlist, mpi_comm_global_h.myid, mpi_comm_global_h.nprocs, true, true);

    for (auto isp = 0; isp != n_spins; isp++)
    {
        Profiler::start("build_real_space_exx_3", "Prepare DM libRI object");
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
        envs::ofs_myid << "Number of Dmat keys: " << get_num_keys(dmat_libri) << "\n";
        // print_keys(envs::ofs_myid, dmat_libri);
        exx_libri.set_Ds(dmat_libri, Params::libri_exx_threshold_D);
        Profiler::stop("build_real_space_exx_3");
        utils::lib_printf("Task %4d: DM setup for EXX\n", mpi_comm_global_h.myid);

        Profiler::start("build_real_space_exx_4", "Call libRI Hexx calculation");
        exx_libri.cal_Hs();
        Profiler::stop("build_real_space_exx_4");

        utils::lib_printf("Task %4d: cal_Hs elapsed time: %f\n", mpi_comm_global_h.myid, Profiler::get_wall_time_last("build_real_space_exx_4"));
        envs::ofs_myid << "Number of exx_libri.Hs keys: " << get_num_keys(exx_libri.Hs) << "\n";
        // print_keys(envs::ofs_myid, exx_libri.Hs);
        // ofs_myid << "exx_libri.Hs:\n" << exx_libri.Hs << endl;

        for (const auto &I_JR_exx: exx_libri.Hs)
        {
            const auto &I = I_JR_exx.first;
            const auto &n_I = atomic_basis_wfc.get_atom_nb(I);
            for (const auto &JR_exx: I_JR_exx.second)
            {
                const auto &J = JR_exx.first.first;
                const auto &n_J = atomic_basis_wfc.get_atom_nb(J);
                const auto &Ra = JR_exx.first.second;
                const auto R = Vector3_Order<int>{Ra[0], Ra[1], Ra[2]};
                Matd exx_temp(n_I, n_J, JR_exx.second.ptr(), MAJOR::ROW);
                this->exx[isp][R][I][J] = exx_temp;
            }
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

#else
    if (mpi_comm_global_h.is_root())
    {
        utils::lib_printf("Error: trying build EXX orbital energy with LibRI, but the program is not compiled against LibRI\n");
    }
    throw std::logic_error("compilation");
    mpi_comm_global_h.barrier();
#endif

    is_rspace_build_= true;
}


void Exx::build_KS(const std::vector<std::vector<ComplexMatrix>> &wfc_target,
                  const std::vector<Vector3_Order<double>> &kfrac_target)
{
    using LIBRPA::envs::mpi_comm_global_h;
    using LIBRPA::envs::blacs_ctxt_global_h;
    using RI::Communicate_Tensors_Map_Judge::comm_map2_first;

    assert(this->is_rspace_build_);
    // Reset k-space matrices built from last call
    if (this->is_kspace_built_)
    {
        utils::lib_printf("Warning: reset EXX k-space matrices\n");
        this->reset_kspace();
    }

    const auto& n_aos = this->mf_.get_n_aos();
    const auto& n_spins = this->mf_.get_n_spins();
    const auto& n_bands = this->mf_.get_n_bands();

    // prepare scalapack array descriptors
    Array_Desc desc_nao_nao(blacs_ctxt_global_h);
    Array_Desc desc_nband_nao(blacs_ctxt_global_h);
    Array_Desc desc_nband_nband(blacs_ctxt_global_h);
    Array_Desc desc_nband_nband_fb(blacs_ctxt_global_h);

    desc_nao_nao.init_1b1p(n_aos, n_aos, 0, 0);
    desc_nband_nao.init_1b1p(n_bands, n_aos, 0, 0);
    desc_nband_nband.init_1b1p(n_bands, n_bands, 0, 0);
    desc_nband_nband_fb.init(n_bands, n_bands, n_bands, n_bands, 0, 0);

    // local 2D-block submatrices
    auto Hexx_nao_nao = init_local_mat<complex<double>>(desc_nao_nao, MAJOR::COL);
    auto temp_nband_nao = init_local_mat<complex<double>>(desc_nband_nao, MAJOR::COL);
    auto Hexx_nband_nband = init_local_mat<complex<double>>(desc_nband_nband, MAJOR::COL);
    auto Hexx_nband_nband_fb = init_local_mat<complex<double>>(desc_nband_nband_fb, MAJOR::COL);

    const auto set_IJ_naonao = LIBRPA::utils::get_necessary_IJ_from_block_2D(
        atomic_basis_wfc, atomic_basis_wfc, desc_nao_nao);
    const auto Iset_Jset = convert_IJset_to_Iset_Jset(set_IJ_naonao);

    for (int isp = 0; isp < n_spins; isp++)
    {
        // collect necessary data
        Profiler::start("build_real_space_exx_5", "Collect Hexx IJ from world");
        const auto &exx_is = this->exx.at(isp);
        std::map<int, std::map<std::pair<int, std::array<int, 3>>, RI::Tensor<double>>> exx_I_JR_local;

        for (const auto &R_IJ_exx: exx_is)
        {
            const auto R = R_IJ_exx.first;
            for (const auto &I_J_exx: R_IJ_exx.second)
            {
                const auto I = I_J_exx.first;
                const auto &n_I = atomic_basis_wfc.get_atom_nb(I);
                for (const auto &J_exx: I_J_exx.second)
                {
                    const auto J = J_exx.first;
                    const auto &n_J = atomic_basis_wfc.get_atom_nb(J);
                    const std::array<int, 3> Ra{R.x, R.y, R.z};
                    exx_I_JR_local[I][{J, Ra}] = RI::Tensor<double>({n_I, n_J}, J_exx.second.sptr());
                }
            }
        }

        // Collect the IJ pair of Hs with all R for Fourier transform
        auto exx_I_JR = comm_map2_first(mpi_comm_global_h.comm,exx_I_JR_local, Iset_Jset.first, Iset_Jset.second);
        exx_I_JR_local.clear();

        // Convert each <I,<J, R>> pair to the nearest neighbour to speed up later Fourier transform
        // while keep the accuracy in further band interpolation.
        // Reuse the cleared-up exx_I_JR_local object
        if (coord_frac.size() > 0)
        {
            for (auto &I_exxJR: exx_I_JR)
            {
                const auto &I = I_exxJR.first;
                for (auto &JR_exx: I_exxJR.second)
                {
                    const auto &J = JR_exx.first.first;
                    const auto &R = JR_exx.first.second;

                    auto distsq = std::numeric_limits<double>::max();
                    Vector3<int> R_IJ;
                    std::array<int, 3> R_bvk;
                    for (int i = -1; i < 2; i++)
                    {
                        R_IJ.x = i * this->period_.x + R[0];
                        for (int j = -1; j < 2; j++)
                        {
                            R_IJ.y = j * this->period_.y + R[1];
                            for (int k = -1; k < 2; k++)
                            {
                                R_IJ.z = k * this->period_.z + R[2];
                                const auto diff =
                                    (Vector3<double>(coord_frac[I][0], coord_frac[I][1],
                                                     coord_frac[I][2]) -
                                     Vector3<double>(coord_frac[J][0], coord_frac[J][1],
                                                     coord_frac[J][2]) -
                                     Vector3<double>(R_IJ.x, R_IJ.y, R_IJ.z)) * latvec;
                                const auto norm2 = diff.norm2();
                                if (norm2 < distsq)
                                {
                                    distsq = norm2;
                                    R_bvk[0] = R_IJ.x;
                                    R_bvk[1] = R_IJ.y;
                                    R_bvk[2] = R_IJ.z;
                                }
                            }
                        }
                    }
                    exx_I_JR_local[I][{J, R_bvk}] = std::move(JR_exx.second);
                }
            }
        }
        else
        {
            exx_I_JR_local = std::move(exx_I_JR);
        }

        exx_I_JR.clear();
        Profiler::stop("build_real_space_exx_5");

        utils::lib_printf("Task %4d: tensor communicate elapsed time: %f\n", mpi_comm_global_h.myid, Profiler::get_wall_time_last("build_real_space_exx_5"));
        // cout << I_JallR_Hs << endl;

        for (int ik = 0; ik < kfrac_target.size(); ik++)
        {
            Hexx_nao_nao.zero_out();
            Profiler::start("build_real_space_exx_6", "Hexx IJ -> 2D block");
            const auto& kfrac = kfrac_target[ik];
            const std::function<complex<double>(const int &, const std::pair<int, std::array<int, 3>> &)>
                fourier = [kfrac](const int &I, const std::pair<int, std::array<int, 3>> &J_Ra)
                {
                    const auto &Ra = J_Ra.second;
                    Vector3<double> R_IJ(Ra[0], Ra[1], Ra[2]);
                    const auto ang = (kfrac * R_IJ) * TWO_PI;
                    return complex<double>{std::cos(ang), std::sin(ang)};
                };
            collect_block_from_IJ_storage_tensor_transform(Hexx_nao_nao, desc_nao_nao, 
                    atomic_basis_wfc, atomic_basis_wfc, fourier, exx_I_JR_local);
            Profiler::stop("build_real_space_exx_6");
            // utils::lib_printf("%s\n", str(Hexx_nao_nao).c_str());
            const auto &wfc_isp_k = wfc_target[isp][ik];
            blacs_ctxt_global_h.barrier();
            const auto wfc_block = get_local_mat(wfc_isp_k.c, MAJOR::ROW, desc_nband_nao, MAJOR::COL).conj();
            // utils::lib_printf("%s\n", str(wfc_block).c_str());
            // utils::lib_printf("%s\n", desc_nao_nao.info_desc().c_str());
            // utils::lib_printf("%s\n", desc_nband_nao.info_desc().c_str());
            Profiler::start("build_real_space_exx_7", "Rotate Hexx ij -> KS");
            ScalapackConnector::pgemm_f('N', 'N', n_bands, n_aos, n_aos, 1.0,
                                        wfc_block.ptr(), 1, 1, desc_nband_nao.desc,
                                        Hexx_nao_nao.ptr(), 1, 1, desc_nao_nao.desc,
                                        0.0,
                                        temp_nband_nao.ptr(), 1, 1, desc_nband_nao.desc);
            ScalapackConnector::pgemm_f('N', 'C', n_bands, n_bands, n_aos, -1.0,
                                        temp_nband_nao.ptr(), 1, 1, desc_nband_nao.desc,
                                        wfc_block.ptr(), 1, 1, desc_nband_nao.desc,
                                        0.0,
                                        Hexx_nband_nband.ptr(), 1, 1, desc_nband_nband.desc);
            Profiler::stop("build_real_space_exx_7");

            // collect to master
            Profiler::start("build_real_space_exx_8", "Collect Eexx to root process");
            ScalapackConnector::pgemr2d_f(n_bands, n_bands,
                                          Hexx_nband_nband.ptr(), 1, 1, desc_nband_nband.desc,
                                          Hexx_nband_nband_fb.ptr(), 1, 1, desc_nband_nband_fb.desc,
                                          desc_nband_nband_fb.ictxt());
            this->exx_is_ik_KS[isp][ik] = Hexx_nband_nband_fb.copy();
            // cout << "Hexx_nband_nband_fb isp " << isp  << " ik " << ik << endl << Hexx_nband_nband_fb;
            if (blacs_ctxt_global_h.myid == 0)
            {
                for (int ib = 0; ib != n_bands; ib++)
                    this->Eexx[isp][ik][ib] = Hexx_nband_nband_fb(ib, ib).real();
            }
            Profiler::stop("build_real_space_exx_8");
        }
    }
}

void Exx::build_KS_kgrid()
{
    this->build_KS(this->mf_.get_eigenvectors(), this->kfrac_list_);
}
void Exx::build_KS_kgrid0()
{
    this->build_KS(this->mf_.get_eigenvectors0(), this->kfrac_list_);
}
void Exx::build_KS_band(const std::vector<std::vector<ComplexMatrix>> &wfc_band,
                        const std::vector<Vector3_Order<double>> &kfrac_band)
{
    this->build_KS(wfc_band, kfrac_band);
}

void Exx::reset_rspace()
{
    this->exx.clear();
    this->is_rspace_build_ = false;
}

void Exx::reset_kspace()
{
    this->exx_is_ik_KS.clear();
    this->Eexx.clear();
    this->is_kspace_built_ = false;
}


}
