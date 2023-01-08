#include "atomic_basis.h"
#include "gw.h"
#include "epsilon.h"
#include "input.h"
#include "libri_utils.h"
#ifdef __USE_LIBRI
#include <RI/physics/GW.h>
#endif

namespace LIBRPA
{

G0W0::G0W0(const MeanField &mf,
             const vector<Vector3_Order<double>>& kfrac_list,
             const TFGrids &tfg)
    : mf(mf), kfrac_list(kfrac_list), tfg(tfg)
{}

void G0W0::build_spacetime(
    const atpair_R_mat_t &LRI_Cs,
    const map<double,
              atom_mapping<std::map<Vector3_Order<double>, matrix_m<complex<double>>>>::pair_t_old>
        &Wc_freq_q,
    const vector<Vector3_Order<int>> &Rlist, const Vector3_Order<int> &R_period)
{
    if (mpi_comm_world_h.myid == 0)
    {
        printf("Calculating correlation self-energy by space-time method\n");
    }
    if (!tfg.has_time_grids())
    {
        if (mpi_comm_world_h.myid == 0)
        {
            printf("Parsed time-frequency object do not have time grids, exiting\n");
        }
        mpi_comm_world_h.barrier();
        throw std::logic_error("no time grids");
    }
#ifndef __USE_LIBRI
    if (mpi_comm_world_h.myid == 0)
    {
        cout << "LIBRA::G0W0::build_spacetime is only implemented on top of LibRI" << endl;
        cout << "Please recompiler libRPA with -DUSE_LIBRI and configure include path" << endl;
    }
    mpi_comm_world_h.barrier();
    throw std::logic_error("compilation");
#else
    const auto Wc_tau_R = CT_FT_Wc_freq_q(Wc_freq_q, tfg, Rlist);
    RI::G0W0<int, int, 3, double> g0w0_libri;
    map<int,std::array<double,3>> atoms_pos;
    for (int i = 0; i != natom; i++)
        atoms_pos.insert(pair<int, std::array<double, 3>>{i, {0, 0, 0}});
    std::array<int,3> period_array{R_period.x,R_period.y,R_period.z};
    g0w0_libri.set_parallel(MPI_COMM_WORLD, atoms_pos, lat_array, period_array);

    // FIXME: duplicate codes to prepare LibRI objects from chi0
    std::vector<std::pair<atom_t, std::pair<atom_t, Vector3_Order<int>>>> IJRs_local;
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
                    IJRs_local.push_back({I, {J, R}});
                    n_Cs_IJRs_local++;
                }
            }
    if (mpi_comm_world_h.myid == 0)
        printf("Total count of Cs: %zu\n", n_Cs_IJRs);
    printf("| Number of Cs on Proc %4d: %zu\n", mpi_comm_world_h.myid, n_Cs_IJRs_local);

    std::map<int, std::map<std::pair<int,std::array<int,3>>,RI::Tensor<double>>> Cs_libri;
    for (auto &IJR: IJRs_local)
    {
        const auto I = IJR.first;
        const auto J = IJR.second.first;
        const auto R = IJR.second.second;
        const std::array<int,3> Ra{R.x,R.y,R.z};
        const auto mat = transpose(*(LRI_Cs.at(I).at(J).at(R)));
        std::valarray<double> mat_array(mat.c, mat.size);
        std::shared_ptr<std::valarray<double>> mat_ptr = std::make_shared<std::valarray<double>>();
        *mat_ptr=mat_array;
        // Tensor<double> Tmat({size_t((*mat).nr),size_t((*mat).nc)},mat_ptr);
        Cs_libri[I][{J, Ra}] = RI::Tensor<double>({atomic_basis_abf.get_atom_nb(I),
                                                   atomic_basis_wfc.get_atom_nb(I),
                                                   atomic_basis_wfc.get_atom_nb(J)}, mat_ptr);
    }
    g0w0_libri.set_Cs(Cs_libri, 0.0);

    auto IJR_local_gf = dispatch_vector_prod(tot_atpair_ordered, Rlist, mpi_comm_world_h.myid, mpi_comm_world_h.nprocs, true, false);
    std::set<Vector3_Order<int>> Rs_local;
    for (const auto &IJR: IJR_local_gf)
    {
        Rs_local.insert(IJR.second);
    }

    for (int ispin = 0; ispin != mf.get_n_spins(); ispin++)
    {
        for (auto itau = 0; itau != tfg.get_n_grids(); itau++)
        {
            const auto tau = tfg.get_time_nodes()[itau];
            // build the Wc LibRI object. Note <JI> has to be converted from <IJ>
            // by W_IJ(R) = W^*_JI(-R)
            const auto &Wc_tau = Wc_tau_R.at(tau);
            std::map<int, std::map<std::pair<int, std::array<int, 3>>, RI::Tensor<double>>> Wc_libri;
            for (const auto &I_JRWc: Wc_tau)
            {
                const auto &I = I_JRWc.first;
                const auto &nabf_I = atomic_basis_abf.get_atom_nb(I);
                for (const auto &J_RWc: I_JRWc.second)
                {
                    const auto &J = J_RWc.first;
                    const auto &nabf_J = atomic_basis_abf.get_atom_nb(J);
                    for (const auto &R_Wc: J_RWc.second)
                    {
                        const auto &R = R_Wc.first;
                        // handle the <IJ(R)> block
                        Wc_libri[I][{J, {R.x, R.y, R.z}}] = RI::Tensor<double>({nabf_I, nabf_J}, R_Wc.second.get_real().sptr());
                        // handle the <JI(R)> block
                        if (I == J) continue;
                        auto minusR = -R % R_period;
                        if (J_RWc.second.count(minusR) == 0) continue;
                        const auto &Wc_IJmR = J_RWc.second.at(minusR).get_real().transpose();
                        Wc_libri[J][{I, {R.x, R.y, R.z}}] = RI::Tensor<double>({nabf_J, nabf_I}, Wc_IJmR.sptr());
                    }
                }
            }

            auto sigc_posi_tau = g0w0_libri.Sigc_tau;
            auto sigc_nega_tau = g0w0_libri.Sigc_tau;

            const auto gf = mf.get_gf_real_imagtimes_Rs(ispin, kfrac_list, {tau, -tau}, {Rs_local.cbegin(), Rs_local.cend()});
            for (auto t: {tau, -tau})
            {
                std::map<int, std::map<std::pair<int, std::array<int, 3>>, RI::Tensor<double>>> gf_libri;
                for (const auto &IJR: IJR_local_gf)
                {
                    const auto &I = IJR.first.first;
                    const auto &n_I = atomic_basis_wfc.get_atom_nb(I);
                    const auto &J = IJR.first.second;
                    const auto &n_J = atomic_basis_wfc.get_atom_nb(J);
                    const auto &R = IJR.second;
                    const auto &gf_global = gf.at(t).at(R);
                    matrix gf_IJ_block(n_I, n_J);
                    {
                        for (int i = 0; i != n_I; i++)
                            for (int j = 0; j != n_J; j++)
                            {
                                gf_IJ_block(i, j) = gf_global(atomic_basis_wfc.get_global_index(I, i),
                                                              atomic_basis_wfc.get_global_index(J, j));
                            }
                        std::shared_ptr<std::valarray<double>> mat_ptr = std::make_shared<std::valarray<double>>(gf_IJ_block.c, gf_IJ_block.size);
                        gf_libri[I][{J, {R.x, R.y, R.z}}] = RI::Tensor<double>({n_I, n_J}, mat_ptr);
                    }
                }
                g0w0_libri.cal_Sigc(gf_libri, 0.0, Wc_libri, 0.0);
                if (t > 0)
                    sigc_posi_tau = std::move(g0w0_libri.Sigc_tau);
                else
                    sigc_nega_tau = std::move(g0w0_libri.Sigc_tau);
            }
            // zmy debug
            // for (const auto &I_JRsigc: sigc_posi_tau)
            // {
            //     for (const auto &JR_sigc: I_JRsigc.second)
            //     {
            //         cout << "I " << I_JRsigc.first << " J(R) " << JR_sigc.first.first << " " << JR_sigc.first.second << endl;
            //         cout << JR_sigc.second;
            //     }
            // }
            // symmetrize and perform transformation
            for (const auto &I_JR_sigc_posi: sigc_posi_tau)
            {
                const auto &I = I_JR_sigc_posi.first;
                const auto &n_I = atomic_basis_wfc.get_atom_nb(I);
                for (const auto &JR_sigc_posi: I_JR_sigc_posi.second)
                {
                    const auto &J = JR_sigc_posi.first.first;
                    const auto &n_J = atomic_basis_wfc.get_atom_nb(J);
                    const auto &Ra = JR_sigc_posi.first.second;
                    const auto &sigc_posi_block = JR_sigc_posi.second;
                    const auto &sigc_nega_block = sigc_nega_tau.at(I).at(JR_sigc_posi.first);
                    auto sigc_cos = sigc_posi_block + sigc_nega_block;
                    auto sigc_sin = sigc_posi_block - sigc_nega_block;
                    for (int ik = 0; ik != mf.get_n_kpoints(); ik++)
                    {
                        const auto &kfrac = kfrac_list[ik];
                        const auto ang = (kfrac * Vector3_Order<double>(Ra[0], Ra[1], Ra[2])) * TWO_PI;
                        const auto kphase = complex<double>{std::cos(ang), std::sin(ang)};
                        for (int iomega = 0; iomega != tfg.get_n_grids(); iomega++)
                        {
                            const auto omega = tfg.get_freq_nodes()[iomega];
                            const auto t2f_sin = tfg.get_sintrans_t2f()(itau, iomega);
                            const auto t2f_cos = tfg.get_costrans_t2f()(itau, iomega);
                            matrix_m<complex<double>> sigc_temp(n_I, n_J, MAJOR::COL);
                            for (int i = 0; i != n_I; i++)
                                for (int j = 0; j != n_J; j++)
                                {
                                    sigc_temp(i, j) += std::complex<double>{sigc_cos(i, j) * t2f_cos, sigc_sin(i, j) * t2f_sin};
                                }
                            sigc_temp *= - kphase;
                            if (sigc_is_freq_k_IJ.count(ispin) == 0 ||
                                sigc_is_freq_k_IJ.at(ispin).count(omega) == 0 ||
                                sigc_is_freq_k_IJ.at(ispin).at(omega).count(kfrac) == 0 ||
                                sigc_is_freq_k_IJ.at(ispin).at(omega).at(kfrac).count(I) == 0 ||
                                sigc_is_freq_k_IJ.at(ispin).at(omega).at(kfrac).at(I).count(J) == 0)
                            {
                                sigc_is_freq_k_IJ[ispin][omega][kfrac][I][J] = std::move(sigc_temp);
                            }
                            else
                            {
                                sigc_is_freq_k_IJ[ispin][omega][kfrac][I][J] += sigc_temp;
                            }
                        }
                    }
                }
            }
        }
    }
#endif
}

} // namespace LIBRPA
