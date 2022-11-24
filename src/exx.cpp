#include "parallel_mpi.h"
#include "params.h"
#include "constants.h"
#include "exx.h"
#include "input.h"
#include "lapack_connector.h"
#include "vector3_order.h"
#ifdef __USE_LIBRI
#include <RI/physics/Exx.h>
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
    if (params.use_libri_exx)
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
#ifdef __USE_LIBRI
    if (mpi_comm_world_h.is_root())
        printf("Computing EXX orbital energy using LibRI\n");
    mpi_comm_world_h.barrier();
    ::Exx<int, int, 3, double> exx_libri;
    map<int,std::array<double,3>> atoms_pos;
    for(int i=0;i!=atom_mu.size();i++)
        atoms_pos.insert(pair<int,std::array<double,3>>{i,{0,0,0}});
    std::array<double,3> xa{latvec.e11,latvec.e12,latvec.e13};
    std::array<double,3> ya{latvec.e21,latvec.e22,latvec.e23};
    std::array<double,3> za{latvec.e31,latvec.e32,latvec.e33};
    std::array<std::array<double,3>,3> lat_array{xa,ya,za};
    std::array<int,3> period_array{R_period.x,R_period.y,R_period.z};
    exx_libri.set_parallel(MPI_COMM_WORLD, atoms_pos, lat_array, period_array);
    exx_libri.set_csm_threshold(params.libri_exx_threshold_CSM);

    // initialize Cs to each process
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
    // if (mpi_comm_world_h.myid == 0)
    //     printf("Total count of Cs: %zu\n", n_IJRs);
    // printf("| Number of Cs on Proc %4d: %zu\n", mpi_comm_world_h.myid, n_IJRs_local);
    std::map<int, std::map<std::pair<int,std::array<int,3>>,Tensor<double>>> Cs_libri;
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
        Cs_libri[I][{J, Ra}] = Tensor<double>({atom_mu[I], atom_nw[I], atom_nw[J]}, mat_ptr);
    }
    exx_libri.set_Cs(Cs_libri, params.libri_exx_threshold_C);

    // initialize Coulomb matrix
    std::map<int, std::map<std::pair<int,std::array<int,3>>,Tensor<double>>> V_libri;
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
        V_libri[I][{J, Ra}] = Tensor<double>({size_t(VIJR->nr), size_t(VIJR->nc)}, pv);
    }
    exx_libri.set_Vs(V_libri, params.libri_exx_threshold_V);

    // initialize density matrix
    vector<atpair_t> atpair_dmat;
    for (int I = 0; I < atom_nw.size(); I++)
        for (int J = 0; J < atom_nw.size(); J++)
            atpair_dmat.push_back({I, J});
    const auto dmat_IJRs_local = dispatch_vector_prod(atpair_dmat, Rlist, mpi_comm_world_h.myid, mpi_comm_world_h.nprocs, true, true);
    for (auto isp = 0; isp != this->mf_.get_n_spins(); isp++)
    {
        std::map<int, std::map<std::pair<int,std::array<int,3>>,Tensor<double>>> dmat_libri;
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
                    dmat_libri[I][{J, Ra}] = Tensor<double>({size_t(dmat_IJR.nr), size_t(dmat_IJR.nc)}, pdmat);
                }
            }
        }
        exx_libri.set_Ds(dmat_libri, params.libri_exx_threshold_D);
        exx_libri.cal_Hs();
        for (auto &T: exx_libri.Hs)
        {
            for (auto &T1: T.second)
            {
                const auto& I = T.first;
                const auto& J = T1.first.first;
                const auto& Ra =  T1.first.second;
                for (const auto& kfrac: this->kfrac_list_)
                {
                    if (Hexx.count(isp) == 0 ||
                        Hexx.at(isp).count(I) == 0 ||
                        Hexx.at(isp).at(I).count(J) == 0 ||
                        Hexx.at(isp).at(I).at(J).count(kfrac) == 0
                        )
                    {
                        Hexx[isp][I][J][kfrac] = make_shared<ComplexMatrix>();
                        ComplexMatrix cm(atom_nw[I], atom_nw[J]);
                        *Hexx[isp][I][J][kfrac] = cm;
                    }
                    ComplexMatrix cm(atom_nw[I], atom_nw[J]);
                    for (int i = 0; i != cm.size; i++)
                        cm.c[i] = *(T1.second.ptr()+i);
                    double ang = kfrac * Vector3_Order<int>{Ra[0], Ra[1], Ra[2]} * TWO_PI;
                    complex<double> kphase = complex<double>(cos(ang), sin(ang));
                    // NOTE: scaling with nkpoints is required
                    *Hexx[isp][I][J][kfrac] += cm * (kphase / double(this->mf_.get_n_kpoints()));
                }
                // debug, check size
                // cout << I << " " << J << " {" << Ra[0] << " " << Ra[1] << " " << Ra[2] << "} whole size: " << T1.second.get_shape_all() << endl;
                // matrix m(atom_nw[I], atom_nw[J]);
                // for (int i = 0; i < m.size; i++)
                //     m.c[i] = (*T1.second.data)[i];
                // print_matrix("", m);
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

    // Rotate the Hamiltonian in basis of KS state for each spin and kpoints
    const auto& n_aos = this->mf_.get_n_aos();
    for (int isp = 0; isp < this->mf_.get_n_spins(); isp++)
    {
        for (int ik = 0; ik < this->mf_.get_n_kpoints(); ik++)
        {
            const auto& k = this->kfrac_list_[ik];
            // retrieve the global Hexx
            ComplexMatrix hexx(n_aos, n_aos);
            for (const auto& I_JkH: this->Hexx.at(isp))
            {
                const auto& I = I_JkH.first;
                const auto I_num = atom_nw[I];
                for (const auto& J_kH: I_JkH.second)
                {
                    const auto& J = J_kH.first;
                    const auto J_num = atom_nw[J];
                    const auto& H = J_kH.second.at(k);
                    for (int i = 0; i != I_num; i++)
                    {
                        size_t i_glo = atom_iw_loc2glo(I, i);
                        for (int j = 0; j != J_num; j++)
                        {
                            size_t j_glo = atom_iw_loc2glo(J, j);
                            hexx(i_glo, j_glo) = (*H)(i, j);
                        }
                    }
                }
            }
            const auto& wfc = this->mf_.get_eigenvectors()[isp][ik];
            this->Hexx_KS[isp][ik] = wfc * hexx * transpose(wfc, true);
            this->Hexx_KS[isp][ik] *= -1.0;
            for (int ib = 0; ib != this->mf_.get_n_bands(); ib++)
                this->Eexx[isp][ik][ib] = this->Hexx_KS[isp][ik](ib, ib).real();
            // debug, print the matrix
            // print_complex_matrix("", this->Hexx_KS[isp][ik] * HA2EV);
        }
    }
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

