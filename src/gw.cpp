#include "gw.h"

#include <fstream>
#include <functional>
#include <map>
#include <sstream>

#include "atomic_basis.h"
#include "constants.h"
#include "envs_blacs.h"
#include "envs_mpi.h"
#include "epsilon.h"
#include "geometry.h"
#include "libri_utils.h"
#include "matrix_m_parallel_utils.h"
#include "params.h"
#include "pbc.h"
#include "profiler.h"
#include "stl_io_helper.h"
#include "utils_blacs.h"
#include "utils_io.h"
#include "utils_mem.h"
#include "utils_mpi_io.h"

#ifdef LIBRPA_USE_LIBRI
#include <RI/global/Tensor.h>
#include <RI/physics/GW.h>
using RI::Tensor;
using RI::Communicate_Tensors_Map_Judge::comm_map2_first;
#endif

namespace LIBRPA
{

G0W0::G0W0(const MeanField &mf, const vector<Vector3_Order<double>> &kfrac_list,
           const TFGrids &tfg_in, const Vector3_Order<int> &period)
    : mf(mf), kfrac_list(kfrac_list), tfg(tfg_in), period_(period)
{
    is_rspace_built_ = false;
    is_kspace_built_ = false;
}

void G0W0::reset_rspace()
{
    sigc_is_f_R_IJ.clear();
    is_rspace_built_ = false;
}

void G0W0::reset_kspace()
{
    sigc_is_ik_f_KS.clear();
    is_kspace_built_ = false;
}

void static unfold_abfs_Wc(
    map<Vector3_Order<double>, ComplexMatrix> &sinvS,
    atom_mapping<std::map<Vector3_Order<double>, matrix_m<complex<double>>>>::pair_t_old &Wc,
    const vector<Vector3_Order<double>> &qlist, map<atom_t, size_t> &atom_mu_large,
    map<atom_t, size_t> &atom_mu_small)
{
    using LIBRPA::Array_Desc;
    using LIBRPA::AtomicBasis;
    using LIBRPA::envs::blacs_ctxt_global_h;
    using LIBRPA::envs::mpi_comm_global_h;

    LIBRPA::atomic_basis_abf.set(atom_mu_small);
    atom_mu_part_range.resize(atom_mu_small.size());
    atom_mu_part_range[0] = 0;
    for (int I = 1; I != atom_mu_small.size(); I++)
        atom_mu_part_range[I] = atom_mu_small.at(I - 1) + atom_mu_part_range[I - 1];
    N_all_mu = atom_mu_part_range[natom - 1] + atom_mu_small[natom - 1];

    // before reset atom_mu: large abfs
    int all_mu = 0;
    vector<int> mu_shift(atom_mu_large.size());
    for (int I = 0; I != atom_mu_large.size(); I++)
    {
        mu_shift[I] = all_mu;
        all_mu += atom_mu_large[I];
    }

    // after reset atom_mu: small abfs
    int all_mu_s = 0;
    vector<int> mu_s_shift(atom_mu_small.size());
    for (int I = 0; I != atom_mu_small.size(); I++)
    {
        mu_s_shift[I] = all_mu_s;
        all_mu_s += atom_mu_small[I];
    }
    AtomicBasis atomic_basis_abf_l;
    atomic_basis_abf_l.set(atom_mu_large);

    const complex<double> CONE{1.0, 0.0};
    Array_Desc desc_nabf_nabf_ll(blacs_ctxt_global_h);
    Array_Desc desc_nabf_nabf_ss(blacs_ctxt_global_h);
    Array_Desc desc_nabf_nabf_sl(blacs_ctxt_global_h);
    desc_nabf_nabf_ll.init_square_blk(all_mu, all_mu, 0, 0);
    desc_nabf_nabf_ss.init_square_blk(all_mu_s, all_mu_s, 0, 0);
    desc_nabf_nabf_sl.init_square_blk(all_mu_s, all_mu, 0, 0);
    const auto set_IJ_nabf_nabf = LIBRPA::utils::get_necessary_IJ_from_block_2D_sy(
        'U', LIBRPA::atomic_basis_abf, desc_nabf_nabf_ss);
    const auto s0_s1 = get_s0_s1_for_comm_map2_first(set_IJ_nabf_nabf);
    auto Wc_block = init_local_mat<complex<double>>(desc_nabf_nabf_ss, MAJOR::COL);
    auto Wcll_block = init_local_mat<complex<double>>(desc_nabf_nabf_ll, MAJOR::COL);
    auto u_block = init_local_mat<complex<double>>(desc_nabf_nabf_sl, MAJOR::COL);
    auto Wc_u = init_local_mat<complex<double>>(desc_nabf_nabf_sl, MAJOR::COL);
    // for 2D->IJ
    int I, iI;
    map<int, vector<int>> map_lor_v;
    map<int, vector<int>> map_loc_v;
    for (int i_lo = 0; i_lo != desc_nabf_nabf_ll.m_loc(); i_lo++)
    {
        int i_glo = desc_nabf_nabf_ll.indx_l2g_r(i_lo);
        atomic_basis_abf_l.get_local_index(i_glo, I, iI);
        map_lor_v[I].push_back(iI);
    }
    for (int i_lo = 0; i_lo != desc_nabf_nabf_ll.n_loc(); i_lo++)
    {
        int i_glo = desc_nabf_nabf_ll.indx_l2g_c(i_lo);
        atomic_basis_abf_l.get_local_index(i_glo, I, iI);
        map_loc_v[I].push_back(iI);
    }

    // IJ pair of shrinked chi0 to be returned
    pair<set<int>, set<int>> Iset_Jset_c;
    const auto atpair_local = dispatch_upper_trangular_tasks(
        natom, blacs_ctxt_global_h.myid, blacs_ctxt_global_h.nprows, blacs_ctxt_global_h.npcols,
        blacs_ctxt_global_h.myprow, blacs_ctxt_global_h.mypcol);
    for (const auto &ap : atpair_local)
    {
        Iset_Jset_c.first.insert(ap.first);
        Iset_Jset_c.second.insert(ap.second);
    }

    for (int iq = 0; iq < qlist.size(); iq++)
    {
        const auto &q = qlist[iq];
        std::array<double, 3> qa = {q.x, q.y, q.z};
        const auto &U = sinvS.at(q);
        // Profiler::start("unfold_prepare_Wc_2d", "Prepare Wc 2D block for unfold");
        Wc_block.zero_out();
        Wcll_block.zero_out();
        u_block.zero_out();
        Wc_u.zero_out();
        {
            std::map<int,
                     std::map<std::pair<int, std::array<double, 3>>, RI::Tensor<complex<double>>>>
                wc_libri;
            atom_mapping<ComplexMatrix>::pair_t_old Wc_IJ;

            for (const auto &IJqc : Wc)
            {
                const auto &I = IJqc.first;
                for (const auto &Jqc : IJqc.second)
                {
                    const auto &J = Jqc.first;
                    if (!Wc_IJ[I].count(J))
                    {
                        Wc_IJ[I][J].create(atom_mu_small[I], atom_mu_small[J]);
                    }
                    for (const auto &qc : Jqc.second)
                    {
                        const auto &qq = qc.first;
                        if (qq == q)
                        {
                            const auto &c = qc.second;
                            for (int ir = 0; ir < c.nr(); ir++)
                            {
                                for (int ic = 0; ic < c.nc(); ic++)
                                {
                                    Wc_IJ[I][J](ir, ic) = c(ir, ic);
                                }
                            }
                        }
                    }
                }
            }

            for (const auto &M_Nchi : Wc_IJ)
            {
                const auto &M = M_Nchi.first;
                const auto n_mu = LIBRPA::atomic_basis_abf.get_atom_nb(M);
                for (const auto &N_chi : M_Nchi.second)
                {
                    const auto &N = N_chi.first;
                    const auto n_nu = LIBRPA::atomic_basis_abf.get_atom_nb(N);
                    const auto &chi = N_chi.second;
                    std::valarray<complex<double>> chi_va(chi.c, chi.size);
                    auto pchi = std::make_shared<std::valarray<complex<double>>>();
                    *pchi = chi_va;
                    wc_libri[M][{N, qa}] = RI::Tensor<complex<double>>({n_mu, n_nu}, pchi);
                }
            }

            // wait for all mpi to calculate chi0_libri
            // then collect chi0_libri to chi0_block
            mpi_comm_global_h.barrier();
            // Profiler::start("unfold_prepare_Wc_2d_comm_map2");
            const auto IJq_wc = RI::Communicate_Tensors_Map_Judge::comm_map2_first(
                mpi_comm_global_h.comm, wc_libri, s0_s1.first, s0_s1.second);
            // Profiler::stop("unfold_prepare_Wc_2d_comm_map2");
            // Profiler::start("unfold_prepare_Wc_2d_collect_block");
            collect_block_from_ALL_IJ_Tensor(Wc_block, desc_nabf_nabf_ss, LIBRPA::atomic_basis_abf,
                                             qa, true, CONE, IJq_wc, MAJOR::ROW);
            // Profiler::stop("unfold_prepare_Wc_2d_collect_block");
        }
        // Profiler::stop("unfold_prepare_Wc_2d");
        for (int ir = 0; ir < U.nr; ir++)
        {
            const int ilo = desc_nabf_nabf_sl.indx_g2l_r(ir);
            if (ilo < 0) continue;
            for (int ic = 0; ic < U.nc; ic++)
            {
                const int jlo = desc_nabf_nabf_sl.indx_g2l_c(ic);
                if (jlo < 0) continue;
                u_block(ilo, jlo) = U(ir, ic);
            }
        }
        // Shape of u_block is N_small x N_large
        ScalapackConnector::pgemm_f('N', 'N', all_mu_s, all_mu, all_mu_s, 1.0, Wc_block.ptr(), 1, 1,
                                    desc_nabf_nabf_ss.desc, u_block.ptr(), 1, 1,
                                    desc_nabf_nabf_sl.desc, 0.0, Wc_u.ptr(), 1, 1,
                                    desc_nabf_nabf_sl.desc);
        ScalapackConnector::pgemm_f('C', 'N', all_mu, all_mu, all_mu_s, 1.0, u_block.ptr(), 1, 1,
                                    desc_nabf_nabf_sl.desc, Wc_u.ptr(), 1, 1,
                                    desc_nabf_nabf_sl.desc, 0.0, Wcll_block.ptr(), 1, 1,
                                    desc_nabf_nabf_ll.desc);

        map<int, map<int, matrix_m<complex<double>>>> chi0s_MNmap;
        map_block_to_IJ_storage_new(chi0s_MNmap, atomic_basis_abf_l, map_lor_v, map_loc_v,
                                    Wcll_block, desc_nabf_nabf_ll, MAJOR::ROW);

        std::map<int, std::map<std::pair<int, std::array<double, 3>>, RI::Tensor<complex<double>>>>
            unfold_Wc_libri;
        for (const auto &M_Nc : chi0s_MNmap)
        {
            const auto &M = M_Nc.first;
            const auto n_mu = atom_mu_large[M];
            for (const auto &N_c : M_Nc.second)
            {
                const auto &N = N_c.first;
                const auto n_nu = atom_mu_large[N];
                const auto &c = N_c.second;
                unfold_Wc_libri[M][{N, qa}] = RI::Tensor<complex<double>>({n_mu, n_nu}, c.sptr());
            }
        }
        const auto IJq_chi = RI::Communicate_Tensors_Map_Judge::comm_map2_first(
            mpi_comm_global_h.comm, unfold_Wc_libri, Iset_Jset_c.first, Iset_Jset_c.second);
        for (auto &IJqc : Wc)
        {
            auto I = IJqc.first;
            for (auto &Jqc : IJqc.second)
            {
                auto J = Jqc.first;
                for (auto &qc : Jqc.second)
                {
                    auto qq = qc.first;
                    if (qq != q) continue;
                    Matz matz_Wc(atom_mu_large[I], atom_mu_large[J]);
                    for (int ir = 0; ir < atom_mu_large[I]; ir++)
                    {
                        for (int ic = 0; ic < atom_mu_large[J]; ic++)
                        {
                            matz_Wc(ir, ic) = IJq_chi.at(I).at({J, qa})(ir, ic);
                        }
                    }
                    qc.second = matz_Wc;
                }
            }
        }
    }

    LIBRPA::atomic_basis_abf.set(atom_mu_large);
    atom_mu_part_range.resize(atom_mu_large.size());
    atom_mu_part_range[0] = 0;
    for (int I = 1; I != atom_mu_large.size(); I++)
        atom_mu_part_range[I] = atom_mu_large.at(I - 1) + atom_mu_part_range[I - 1];
    N_all_mu = atom_mu_part_range[natom - 1] + atom_mu_large[natom - 1];
}

template <typename Tdata>
void G0W0::build_spacetime(
    const Cs_LRI &LRI_Cs,
    map<double,
        atom_mapping<std::map<Vector3_Order<double>, matrix_m<complex<double>>>>::pair_t_old>
        &Wc_freq_q,
    const vector<Vector3_Order<int>> &Rlist, const vector<Vector3_Order<double>> &qlist,
    std::map<Vector3_Order<double>, ComplexMatrix> &sinvS)
{
    using LIBRPA::envs::mpi_comm_global_h;

    if (parallel_routing != ParallelRouting::LIBRI)
    {
        mpi_comm_global_h.barrier();
        throw std::logic_error("not implemented");
    }

    LIBRPA::utils::lib_printf_root("Calculating correlation self-energy by space-time method\n");
    if (!tfg.has_time_grids())
    {
        LIBRPA::utils::lib_printf_root(
            "Parsed time-frequency object do not have time grids, exiting\n");
        mpi_comm_global_h.barrier();
        throw std::logic_error("no time grids");
    }
    mpi_comm_global_h.barrier();

    std::ofstream ofs_sigmac_r;

#ifndef LIBRPA_USE_LIBRI
    if (mpi_comm_global_h.myid == 0)
    {
        cout << "LIBRA::G0W0::build_spacetime is only implemented on top of LibRI" << endl;
        cout << "Please recompiler LibRPA with -DUSE_LIBRI and configure include path" << endl;
    }
    mpi_comm_global_h.barrier();
    throw std::logic_error("compilation");
#else
    map<double, atom_mapping<std::map<Vector3_Order<int>, matrix_m<complex<double>>>>::pair_t_old>
        Wc_tau_R;
    map<double,
        atom_mapping<std::map<Vector3_Order<double>, matrix_m<complex<double>>>>::pair_t_old>
        Wc_tau_q;
    // Transform
    if (Params::use_shrink_abfs)
    {
        if (Params::output_Wc_Rf_mat == 1)
        {
            Profiler::start("unfold_Wc_q", "Unfold Wc (q,w=0)");
            const double freq = tfg.get_freq_nodes()[0];
            auto Wc_q_f0 = Wc_freq_q.at(freq);
            unfold_abfs_Wc(sinvS, Wc_q_f0, qlist, atom_mu_l, atom_mu_s);
            Profiler::stop("unfold_Wc_q");
            Profiler::start("construct_Wc_lower_half", "Construct Lower Half of Wc(q,w)");
            // NOTE: only upper half of Wc is built now
            //       here we recover the other half before transform to R space using the Hermitian property
            //       and Hermitize the diagonal blocks (due to numerical noise)
            vector<atom_t> iatoms_row;
            for (const auto &Mu_NuqWc : Wc_q_f0) iatoms_row.push_back(Mu_NuqWc.first);
            for (auto iatom_row: iatoms_row)
            {
                vector<atom_t> iatoms_col;
                for (const auto &Nu_qWc : Wc_q_f0.at(iatom_row))
                {
                    iatoms_col.push_back(Nu_qWc.first);
                }
                for (auto iatom_col : iatoms_col)
                {
                    for (const auto &q_Wc : Wc_q_f0.at(iatom_row).at(iatom_col))
                    {
                        assert(q_Wc.second.major() == MAJOR::ROW);
                        if(iatom_row != iatom_col)
                            Wc_q_f0[iatom_col][iatom_row][q_Wc.first] = q_Wc.second.get_transpose(true);
                        else // Hermitize the diagonal blocks
                        {
                            auto Wc_mat = q_Wc.second;
                            Wc_mat = (Wc_mat + Wc_mat.get_transpose(true)) * 0.5;
                            Wc_q_f0[iatom_row][iatom_row][q_Wc.first] = Wc_mat;
                        }
                    }
                }
            }

            mpi_comm_global_h.barrier();
            Profiler::stop("construct_Wc_lower_half");
            FT_Wc_q2R(Wc_q_f0, tfg, meanfield.get_n_kpoints(), Rlist, true);
            Wc_q_f0.clear();
        }
        else if (Params::output_Wc_Rf_mat > 1)
        {
            throw std::logic_error("use_shrink_abfs doesn't support output_Wc_Rf_mat > 1");
        }
        Profiler::start("g0w0_build_spacetime_wt_ft_wc", "Tranform Wc (q,w) -> (q,t)");
        Wc_tau_q = CT_Wc_freq2time_q(Wc_freq_q, tfg, meanfield.get_n_kpoints(), Rlist, qlist);

        // HACK: Free up Wc_freq_q to save memory, especially for large Coulomb matrix case and many
        // minimax grids
        Wc_freq_q.clear();
        utils::release_free_mem();
        Profiler::stop("g0w0_build_spacetime_wt_ft_wc");
        LIBRPA::utils::lib_printf_root(
            "Time for Fourier transform of Wc in GW (seconds, Wall/CPU): %f %f\n",
            Profiler::get_wall_time_last("g0w0_build_spacetime_wt_ft_wc"),
            Profiler::get_cpu_time_last("g0w0_build_spacetime_wt_ft_wc"));
    }
    else
    {
        Profiler::start("g0w0_build_spacetime_ct_ft_wc", "Tranform Wc (q,w) -> (R,t)");
        if (Params::output_Wc_Rf_mat > 0)
        {
            Wc_tau_R = CT_FT_Wc_q2R_freq2time(Wc_freq_q, tfg, meanfield.get_n_kpoints(), Rlist);
        }
        else
        {
            Wc_tau_R = CT_FT_Wc_freq_q(Wc_freq_q, tfg, meanfield.get_n_kpoints(), Rlist);
            // HACK: Free up Wc_freq_q to save memory, especially for large Coulomb matrix case and many
            // minimax grids
            Wc_freq_q.clear();
            utils::release_free_mem();
        }
        Profiler::stop("g0w0_build_spacetime_ct_ft_wc");
        LIBRPA::utils::lib_printf_root(
            "Time for Fourier transform of Wc in GW (seconds, Wall/CPU): %f %f\n",
            Profiler::get_wall_time_last("g0w0_build_spacetime_ct_ft_wc"),
            Profiler::get_cpu_time_last("g0w0_build_spacetime_ct_ft_wc"));
    }

    RI::GW<int, int, 3, Tdata> gw_libri;
    map<int, std::array<double, 3>> atoms_pos;
    for (int i = 0; i != natom; i++)
        atoms_pos.insert(pair<int, std::array<double, 3>>{i, {0, 0, 0}});
    std::array<int, 3> period_array{this->period_.x, this->period_.y, this->period_.z};

    Profiler::start("g0w0_build_spacetime_2", "Setup LibRI G0W0 object and C data");
    gw_libri.set_parallel(mpi_comm_global_h.comm, atoms_pos, lat_array, period_array);
    // TODO: template Cs_LRI
    if constexpr (std::is_same<Tdata, std::complex<double>>::value)
    {
        std::map<int, std::map<libri_types<int, int>::TAC, RI::Tensor<Tdata>>> data_libri;
        for (const auto &I_JR_C : LRI_Cs.data_libri)
        {
            const auto I = I_JR_C.first;
            for (const auto &JR_C : I_JR_C.second)
            {
                const auto J = JR_C.first.first;
                const auto R = JR_C.first.second;
                const auto &C = JR_C.second;
                auto JR = std::pair<int, std::array<int, 3>>(J, R);
                data_libri[I][JR] = RI::Global_Func::convert<Tdata>(C);
            }
        }
        gw_libri.set_Cs(data_libri, Params::libri_g0w0_threshold_C);
    }
    else
        gw_libri.set_Cs(LRI_Cs.data_libri, Params::libri_g0w0_threshold_C);

    Profiler::stop("g0w0_build_spacetime_2");
    LIBRPA::utils::lib_printf_root("Time for LibRI G0W0 setup (seconds, Wall/CPU): %f %f\n",
                                   profiler.get_wall_time_last("g0w0_build_spacetime_2"),
                                   profiler.get_cpu_time_last("g0w0_build_spacetime_2"));

    auto IJR_local_gf = dispatch_vector_prod(tot_atpair_ordered, Rlist, mpi_comm_global_h.myid,
                                             mpi_comm_global_h.nprocs, true, false);
    envs::ofs_myid << "IJR_local_gf: " << IJR_local_gf << endl;
    std::set<Vector3_Order<int>> Rs_local;
    for (const auto &IJR : IJR_local_gf)
    {
        Rs_local.insert(IJR.second);
    }

    for (auto itau = 0; itau != tfg.get_n_grids(); itau++)
    {
        const auto tau = tfg.get_time_nodes()[itau];
        if (Params::use_shrink_abfs)
        {
            Profiler::start("unfold_Wc_abfs", "Do shrink transformation");
            unfold_abfs_Wc(sinvS, Wc_tau_q[tau], qlist, atom_mu_l, atom_mu_s);
            Profiler::stop("unfold_Wc_abfs");
            Profiler::start("g0w0_build_spacetime_Rq_ft_wc", "Tranform Wc (q,t) -> (R,t)");
            Wc_tau_R[tau] =
                FT_Wc_q2R(Wc_tau_q[tau], tfg, meanfield.get_n_kpoints(), Rlist, false);
            Profiler::stop("g0w0_build_spacetime_Rq_ft_wc");
            Wc_tau_q[tau].clear();
            atom_mu = atom_mu_l;
            LIBRPA::atomic_basis_abf.set(atom_mu);
            atom_mu_part_range.resize(atom_mu.size());
            atom_mu_part_range[0] = 0;
            for (int I = 1; I != atom_mu.size(); I++)
                atom_mu_part_range[I] = atom_mu.at(I - 1) + atom_mu_part_range[I - 1];

            N_all_mu = atom_mu_part_range[natom - 1] + atom_mu[natom - 1];
        }
        Profiler::start("g0w0_build_spacetime_3", "Prepare LibRI Wc object");
        // LIBRPA::utils::lib_printf("task %d itau %d start\n", mpi_comm_global_h.myid, itau);
        // build the Wc LibRI object. Note <JI> has to be converted from <IJ>
        // by W_IJ(R) = W^*_JI(-R)
        std::map<int, std::map<std::pair<int, std::array<int, 3>>, RI::Tensor<Tdata>>> Wc_libri;
        // in R-tau routing, some process can have zero time point
        // check to avoid out-of-range by at()
        // LIBRPA::utils::lib_printf("task %d Wc_tau_R.count(tau) %zu\n", mpi_comm_global_h.myid,
        // Wc_tau_R.count(tau));
        if (Wc_tau_R.count(tau))
        {
            for (const auto &I_JRWc : Wc_tau_R.at(tau))
            {
                const auto &I = I_JRWc.first;
                const auto &nabf_I = atomic_basis_abf.get_atom_nb(I);
                for (const auto &J_RWc : I_JRWc.second)
                {
                    const auto &J = J_RWc.first;
                    const auto &nabf_J = atomic_basis_abf.get_atom_nb(J);
                    for (const auto &R_Wc : J_RWc.second)
                    {
                        const auto &R = R_Wc.first;
                        // handle the <IJ(R)> block
                        if constexpr (std::is_same<Tdata, std::complex<double>>::value)
                            Wc_libri[static_cast<int>(I)][{static_cast<int>(J), {R.x, R.y, R.z}}] =
                                RI::Tensor<std::complex<double>>({nabf_I, nabf_J},
                                                                 R_Wc.second.sptr());
                        else
                            Wc_libri[static_cast<int>(I)][{static_cast<int>(J), {R.x, R.y, R.z}}] =
                                RI::Tensor<double>({nabf_I, nabf_J}, R_Wc.second.get_real().sptr());
                        // cout << "I " << I << " J " << J <<  " R " << R << " tau " << tau << endl; 
                        // cout << Wc_libri[static_cast<int>(I)][{static_cast<int>(J), {R.x, R.y, R.z}}] << endl; 
                        // std::cout << "R_Wc.second: " << R_Wc.second << std::endl;
                        // handle the <JI(R)> block
                        if (Params::output_Wc_Rf_mat > 0 && !Params::use_shrink_abfs) continue; // full atom-pair has been constructed
                        if (I == J) continue;
                        auto minusR = (-R) % this->period_;
                        if (J_RWc.second.count(minusR) == 0) continue;
                        if constexpr (std::is_same<Tdata, std::complex<double>>::value)
                        {
                            const auto Wc_IJmR = J_RWc.second.at(minusR).get_transpose(true);
                            Wc_libri[static_cast<int>(J)][{static_cast<int>(I), {R.x, R.y, R.z}}] =
                                RI::Tensor<std::complex<double>>({nabf_J, nabf_I}, Wc_IJmR.sptr());
                        }
                        else
                        {
                            const auto Wc_IJmR = J_RWc.second.at(minusR).get_real().get_transpose();
                            Wc_libri[static_cast<int>(J)][{static_cast<int>(I), {R.x, R.y, R.z}}] =
                                RI::Tensor<double>({nabf_J, nabf_I}, Wc_IJmR.sptr());
                        }
                    }
                }
            }
        }
        size_t n_obj_wc_libri = 0;
        for (const auto &w : Wc_libri)
        {
            n_obj_wc_libri += w.second.size();
        }
        profiler.stop("g0w0_build_spacetime_3");
        LIBRPA::utils::lib_printf_root(
            "Time for preparing LibRI Wc object, i_tau %d (seconds, Wall/CPU): %f %f\n", itau + 1,
            profiler.get_wall_time_last("g0w0_build_spacetime_3"),
            profiler.get_cpu_time_last("g0w0_build_spacetime_3"));

        for (int ispin = 0; ispin != mf.get_n_spins(); ispin++)
        {
            for (int isoc1 = 0; isoc1 != mf.get_n_soc(); isoc1++)
            {
                for (int isoc2 = 0; isoc2 != mf.get_n_soc(); isoc2++)
                {
                    std::map<int, std::map<std::pair<int, std::array<int, 3>>, Tensor<Tdata>>>
                        sigc_posi_tau;
                    std::map<int, std::map<std::pair<int, std::array<int, 3>>, Tensor<Tdata>>>
                        sigc_nega_tau;

                    Profiler::start("g0w0_build_spacetime_4", "Compute G(R,t) and G(R,-t)");
                    // NOTE: ``if constexpr`` needs C++-17
                    if constexpr (std::is_same<Tdata, std::complex<double>>::value)
                    {
                        auto gf = mf.get_gf_cplx_imagtimes_Rs(ispin, isoc1, isoc2, kfrac_list,
                                                              {tau, -tau},
                                                              {Rs_local.cbegin(), Rs_local.cend()});
                        std::map<double, std::map<int, std::map<std::pair<int, std::array<int, 3>>,
                                                                RI::Tensor<Tdata>>>>
                            tau_gf_libri;
                        for (auto t : {tau, -tau})
                        {
                            tau_gf_libri[t] = {};
                        }
                        for (const auto &IJR : IJR_local_gf)
                        {
                            const auto &I = IJR.first.first;
                            const auto &n_I = atomic_basis_wfc.get_atom_nb(I);
                            const auto &J = IJR.first.second;
                            const auto &n_J = atomic_basis_wfc.get_atom_nb(J);
                            const auto &R = IJR.second;
                            for (auto t : {tau, -tau})
                            {
                                const auto &gf_global = gf.at(t).at(R);
                                ComplexMatrix gf_IJ_block(n_I, n_J);
                                {
                                    for (int i = 0; i != n_I; i++)
                                        for (int j = 0; j != n_J; j++)
                                        {
                                            gf_IJ_block(i, j) =
                                                gf_global(atomic_basis_wfc.get_global_index(I, i),
                                                          atomic_basis_wfc.get_global_index(J, j));
                                        }
                                    std::shared_ptr<std::valarray<Tdata>> mat_ptr =
                                        std::make_shared<std::valarray<Tdata>>(gf_IJ_block.c,
                                                                               gf_IJ_block.size);
                                    tau_gf_libri[t][static_cast<int>(I)]
                                                [{static_cast<int>(J), {R.x, R.y, R.z}}] =
                                                    RI::Tensor<Tdata>({n_I, n_J}, mat_ptr);
                                }
                            }
                        }
                        gf.clear();
                        Profiler::stop("g0w0_build_spacetime_4");
                        LIBRPA::utils::lib_printf_root(
                            "Time for Green's function, i_spin %d i_tau %d (seconds, Wall/CPU): %f "
                            "%f\n",
                            ispin + 1, itau + 1,
                            Profiler::get_wall_time_last("g0w0_build_spacetime_4"),
                            Profiler::get_cpu_time_last("g0w0_build_spacetime_4"));

                        for (auto t : {tau, -tau})
                        {
                            const auto &gf_libri = tau_gf_libri.at(t);
                            size_t n_obj_gf_libri = 0;
                            for (const auto &gf : gf_libri) n_obj_gf_libri += gf.second.size();

                            double wtime_g0w0_cal_sigc = omp_get_wtime();
                            gw_libri.set_Gs(gf_libri, Params::libri_g0w0_threshold_G);
                            Profiler::start("g0w0_build_spacetime_5", "Call libRI cal_Sigc");
                            gw_libri.cal_Sigmas();
                            Profiler::stop("g0w0_build_spacetime_5");
                            Profiler::start("g0w0_build_spacetime_5_clean");
                            gw_libri.free_Gs();
                            Profiler::stop("g0w0_build_spacetime_5_clean");

                            // Check size of data
                            double mem_mb = get_tensor_map_bytes(gw_libri.Sigmas) * 1e-6;
                            envs::ofs_myid << "Temporary Sigc_tau size for time " << t
                                           << " [MB]: " << mem_mb << endl;

                            if (t > 0)
                                sigc_posi_tau = std::move(gw_libri.Sigmas);
                            else
                                sigc_nega_tau = std::move(gw_libri.Sigmas);
                            gw_libri.Sigmas.clear();

                            wtime_g0w0_cal_sigc = omp_get_wtime() - wtime_g0w0_cal_sigc;
                            LIBRPA::utils::lib_printf(
                                "Task %4d. libRI G0W0, spin %1d, time grid %12.6f. Wc size %zu, GF "
                                "size "
                                "%zu. "
                                "Wall time %f\n",
                                mpi_comm_global_h.myid, ispin, t, n_obj_wc_libri, n_obj_gf_libri,
                                wtime_g0w0_cal_sigc);
                        }
                    }
                    else  // Non-SOC
                    {
                        auto gf = mf.get_gf_real_imagtimes_Rs(ispin, isoc1, isoc2, kfrac_list,
                                                              {tau, -tau},
                                                              {Rs_local.cbegin(), Rs_local.cend()});
                        std::map<double, std::map<int, std::map<std::pair<int, std::array<int, 3>>,
                                                                RI::Tensor<Tdata>>>>
                            tau_gf_libri;
                        for (auto t : {tau, -tau})
                        {
                            tau_gf_libri[t] = {};
                        }
                        for (const auto &IJR : IJR_local_gf)
                        {
                            const auto &I = IJR.first.first;
                            const auto &n_I = atomic_basis_wfc.get_atom_nb(I);
                            const auto &J = IJR.first.second;
                            const auto &n_J = atomic_basis_wfc.get_atom_nb(J);
                            const auto &R = IJR.second;
                            for (auto t : {tau, -tau})
                            {
                                if (gf.count(t) && gf.at(t).count(R))
                                {
                                    const auto &gf_global = gf.at(t).at(R);
                                    matrix gf_IJ_block(n_I, n_J);
                                    {
                                        for (int i = 0; i != n_I; i++)
                                            for (int j = 0; j != n_J; j++)
                                            {
                                                gf_IJ_block(i, j) = gf_global(
                                                    atomic_basis_wfc.get_global_index(I, i),
                                                    atomic_basis_wfc.get_global_index(J, j));
                                            }
                                        std::shared_ptr<std::valarray<Tdata>> mat_ptr =
                                            std::make_shared<std::valarray<Tdata>>(
                                                gf_IJ_block.c, gf_IJ_block.size);

                                        tau_gf_libri[t][static_cast<int>(I)]
                                                    [{static_cast<int>(J), {R.x, R.y, R.z}}] =
                                                        RI::Tensor<Tdata>({n_I, n_J}, mat_ptr);
                                    }
                                }
                            }
                        }
                        gf.clear();
                        Profiler::stop("g0w0_build_spacetime_4");
                        LIBRPA::utils::lib_printf_root(
                            "Time for Green's function, i_spin %d i_tau %d (seconds, Wall/CPU): %f "
                            "%f\n",
                            ispin + 1, itau + 1,
                            Profiler::get_wall_time_last("g0w0_build_spacetime_4"),
                            Profiler::get_cpu_time_last("g0w0_build_spacetime_4"));

                        for (auto t : {tau, -tau})
                        {
                            const auto &gf_libri = tau_gf_libri.at(t);
                            size_t n_obj_gf_libri = 0;
                            for (const auto &gf : gf_libri) n_obj_gf_libri += gf.second.size();

                            double wtime_g0w0_cal_sigc = omp_get_wtime();
                            gw_libri.set_Gs(gf_libri, Params::libri_g0w0_threshold_G);
                            Profiler::start("g0w0_build_spacetime_5", "Call libRI cal_Sigc");
                            gw_libri.cal_Sigmas();
                            Profiler::stop("g0w0_build_spacetime_5");
                            Profiler::start("g0w0_build_spacetime_5_clean");
                            gw_libri.free_Gs();
                            Profiler::stop("g0w0_build_spacetime_5_clean");

                            // Check size of data
                            double mem_mb = get_tensor_map_bytes(gw_libri.Sigmas) * 1e-6;
                            envs::ofs_myid << "Temporary Sigc_tau size for time " << t
                                           << " [MB]: " << mem_mb << endl;

                            if (t > 0)
                                sigc_posi_tau = std::move(gw_libri.Sigmas);
                            else
                                sigc_nega_tau = std::move(gw_libri.Sigmas);
                            gw_libri.Sigmas.clear();

                            wtime_g0w0_cal_sigc = omp_get_wtime() - wtime_g0w0_cal_sigc;
                            LIBRPA::utils::lib_printf(
                                "Task %4d. libRI G0W0, spin %1d, time grid %12.6f. Wc size %zu, GF "
                                "size "
                                "%zu. "
                                "Wall time %f\n",
                                mpi_comm_global_h.myid, ispin, t, n_obj_wc_libri, n_obj_gf_libri,
                                wtime_g0w0_cal_sigc);
                        }
                    }

                    size_t n_IJR_myid = 0;  // for sigcmat output

                    if (Params::output_gw_sigc_mat_rt)
                    {
                        std::stringstream ss;
                        ss << Params::output_dir << "SigcRT"
                           << "_ispin_" << std::setfill('0') << std::setw(2) << ispin << "_s_"
                           << std::setw(1) << isoc1 << std::setw(1) << isoc2 << "_itau_"
                           << std::setfill('0') << std::setw(3) << itau << "_myid_"
                           << std::setfill('0') << std::setw(5) << envs::myid_global << ".dat";
                        ofs_sigmac_r.open(ss.str(), std::ios::out | std::ios::binary);
                        ofs_sigmac_r.write((char *)&n_IJR_myid, sizeof(size_t));  // placeholder
                    }

                    // symmetrize and perform transformation
                    Profiler::start("g0w0_build_spacetime_6", "Transform Sigc (R,t) -> (R,w)");
                    for (const auto &I_JR_sigc_posi : sigc_posi_tau)
                    {
                        const auto &I = I_JR_sigc_posi.first;
                        const auto &n_I = atomic_basis_wfc.get_atom_nb(I);
                        for (const auto &JR_sigc_posi : I_JR_sigc_posi.second)
                        {
                            const auto &J = JR_sigc_posi.first.first;
                            const auto &n_J = atomic_basis_wfc.get_atom_nb(J);
                            const auto &Ra = JR_sigc_posi.first.second;
                            const auto &sigc_posi_block = JR_sigc_posi.second;
                            if (sigc_nega_tau.count(I) == 0 ||
                                sigc_nega_tau.at(I).count(JR_sigc_posi.first) == 0)
                            {
                                continue;
                            }
                            const auto &sigc_nega_block =
                                sigc_nega_tau.at(I).at(JR_sigc_posi.first);
                            const auto R = Vector3_Order<int>(Ra[0], Ra[1], Ra[2]);
                            const auto iR = std::distance(
                                Rlist.cbegin(), std::find(Rlist.cbegin(), Rlist.cend(), R));

                            const auto sigc_cos = Tdata(0.5) * (sigc_posi_block + sigc_nega_block);
                            const auto sigc_sin = Tdata(0.5) * (sigc_posi_block - sigc_nega_block);

                            n_IJR_myid++;
                            if (Params::output_gw_sigc_mat_rt)
                            {
                                // write indices and dimensions
                                size_t dims[5];
                                dims[0] = iR;
                                dims[1] = I;
                                dims[2] = J;
                                dims[3] = n_I;
                                dims[4] = n_J;
                                ofs_sigmac_r.write((char *)dims, 5 * sizeof(size_t));
                                // cos: contribute to the real part
                                ofs_sigmac_r.write((char *)sigc_cos.ptr(),
                                                   n_I * n_J * sizeof(double));
                                // sin: contribute to the imaginary part
                                ofs_sigmac_r.write((char *)sigc_sin.ptr(),
                                                   n_I * n_J * sizeof(double));
                            }

                            for (int iomega = 0; iomega != tfg.get_n_grids(); iomega++)
                            {
                                const auto omega = tfg.get_freq_nodes()[iomega];
                                const auto t2f_sin = tfg.get_sintrans_t2f()(iomega, itau);
                                const auto t2f_cos = tfg.get_costrans_t2f()(iomega, itau);
                                // row-major used here for libRI communication when building
                                // sigc_KS
                                matrix_m<complex<double>> sigc_temp(n_I, n_J, MAJOR::ROW);
                                matrix_m<complex<double>> sigc_temp_minus(n_I, n_J, MAJOR::ROW);
                                for (int i = 0; i != n_I; i++)
                                {
                                    if constexpr (std::is_same<Tdata, std::complex<double>>::value)
                                    {
                                        for (int j = 0; j != n_J; j++)
                                        {
                                            sigc_temp(i, j) = sigc_cos(i, j) * t2f_cos +
                                                              sigc_sin(i, j) * t2f_sin * 1.0i;
                                        }
                                    }
                                    else
                                    {
                                        for (int j = 0; j != n_J; j++)
                                        {
                                            sigc_temp(i, j) = std::complex<double>{
                                                sigc_cos(i, j) * t2f_cos, sigc_sin(i, j) * t2f_sin};
                                            sigc_temp_minus(i, j) = std::complex<double>{sigc_cos(i, j) * t2f_cos, -sigc_sin(i, j) * t2f_sin};
                                        }
                                    }
                                }
                                // TODO: incorporating following code for SCRPA
                                // if (sigc_is_f_R_IJ.count(ispin) == 0 ||
                                //     sigc_is_f_R_IJ.at(ispin).count(-omega) == 0 ||
                                //     sigc_is_f_R_IJ.at(ispin).at(-omega).count(R) == 0 ||
                                //     sigc_is_f_R_IJ.at(ispin).at(-omega).at(R).count(I) == 0 ||
                                //     sigc_is_f_R_IJ.at(ispin).at(-omega).at(R).at(I).count(J) == 0)
                                // {
                                //     sigc_is_f_R_IJ[ispin][-omega][R][I][J] = std::move(sigc_temp_minus);
                                // }
                                // else
                                // {
                                //     sigc_is_f_R_IJ[ispin][-omega][R][I][J] += sigc_temp_minus;
                                // }
                                if (sigc_is_f_R_IJ.count(ispin) == 0 ||
                                    sigc_is_f_R_IJ.at(ispin).count(isoc1) == 0 ||
                                    sigc_is_f_R_IJ.at(ispin).at(isoc1).count(isoc2) == 0 ||
                                    sigc_is_f_R_IJ.at(ispin).at(isoc1).at(isoc2).count(omega) ==
                                        0 ||
                                    sigc_is_f_R_IJ.at(ispin).at(isoc1).at(isoc2).at(omega).count(
                                        R) == 0 ||
                                    sigc_is_f_R_IJ.at(ispin)
                                            .at(isoc1)
                                            .at(isoc2)
                                            .at(omega)
                                            .at(R)
                                            .count(I) == 0 ||
                                    sigc_is_f_R_IJ.at(ispin)
                                            .at(isoc1)
                                            .at(isoc2)
                                            .at(omega)
                                            .at(R)
                                            .at(I)
                                            .count(J) == 0)
                                {
                                    sigc_is_f_R_IJ[ispin][isoc1][isoc2][omega][R][I][J] =
                                        std::move(sigc_temp);
                                }
                                else
                                {
                                    sigc_is_f_R_IJ[ispin][isoc1][isoc2][omega][R][I][J] +=
                                        sigc_temp;
                                }
                            }
                        }
                    }

                    sigc_posi_tau.clear();
                    sigc_nega_tau.clear();

                    Profiler::stop("g0w0_build_spacetime_6");

                    if (Params::output_gw_sigc_mat_rt)
                    {
                        ofs_sigmac_r.seekp(0);
                        ofs_sigmac_r.write((char *)&n_IJR_myid, sizeof(size_t));  // placeholder
                        ofs_sigmac_r.close();
                    }
                }  // isoc2
            }  // isoc1
        }  // ispin
        Profiler::start("g0w0_build_spacetime_free_Ws");
        gw_libri.free_Ws();
        Profiler::stop("g0w0_build_spacetime_free_Ws");
    }
    is_rspace_built_ = true;
    if (Params::use_shrink_abfs)
    {
        sinvS.clear();
    }
#endif

    // Export real-space imaginary-frequency NAO sigma_c matrices
    if (is_rspace_built_ && Params::output_gw_sigc_mat_rf)
    {
        for (int ispin = 0; ispin != mf.get_n_spins(); ispin++)
        {
            for (int isoc1 = 0; isoc1 != mf.get_n_soc(); isoc1++)
            {
                for (int isoc2 = 0; isoc2 != mf.get_n_soc(); isoc2++)
                {
                    for (auto iomega = 0; iomega != tfg.get_n_grids(); iomega++)
                    {
                        size_t n_IJR_myid = 0;
                        std::stringstream ss;
                        ss << Params::output_dir << "SigcRF"
                           << "_ispin_" << std::setfill('0') << std::setw(2) << ispin << "_s_"
                           << std::setw(1) << isoc1 << std::setw(1) << isoc2 << "_iomega_"
                           << std::setfill('0') << std::setw(3) << iomega << "_myid_"
                           << std::setfill('0') << std::setw(5) << envs::myid_global << ".dat";
                        ofs_sigmac_r.open(ss.str(), std::ios::out | std::ios::binary);
                        ofs_sigmac_r.write((char *)&n_IJR_myid, sizeof(size_t));  // placeholder

                        const auto omega = tfg.get_freq_nodes()[iomega];
                        const auto &sigc_RIJ = sigc_is_f_R_IJ[ispin][isoc1][isoc2][omega];
                        for (const auto &R_IJsigc : sigc_RIJ)
                        {
                            const auto &R = R_IJsigc.first;
                            const auto iR = std::distance(
                                Rlist.cbegin(), std::find(Rlist.cbegin(), Rlist.cend(), R));
                            for (const auto &I_Jsigc : R_IJsigc.second)
                            {
                                const auto &I = I_Jsigc.first;
                                const auto &n_I = atomic_basis_wfc.get_atom_nb(I);
                                for (const auto &J_sigc : I_Jsigc.second)
                                {
                                    const auto &J = J_sigc.first;
                                    const auto &n_J = atomic_basis_wfc.get_atom_nb(J);
                                    const auto &sigc = J_sigc.second;
                                    n_IJR_myid++;
                                    size_t dims[5];
                                    dims[0] = iR;
                                    dims[1] = I;
                                    dims[2] = J;
                                    dims[3] = n_I;
                                    dims[4] = n_J;
                                    assert(sigc.size() == n_I * n_J);
                                    ofs_sigmac_r.write((char *)dims, 5 * sizeof(size_t));
                                    ofs_sigmac_r.write((char *)sigc.ptr(),
                                                       sigc.size() * 2 * sizeof(double));
                                }
                            }
                        }
                        ofs_sigmac_r.seekp(0);
                        ofs_sigmac_r.write((char *)&n_IJR_myid, sizeof(size_t));  // placeholder
                        ofs_sigmac_r.close();
                    }
                }
            }
        }
    }
}

void G0W0::read_sigc(const std::string &input_dir, const vector<Vector3_Order<int>> &Rlist)
{
    using LIBRPA::envs::mpi_comm_global_h;

    LIBRPA::utils::lib_printf_root("Reading real-space imaginary-frequency NAO sigma_c matrices\n");
    Profiler::start("g0w0_read_sigc(R,iw) in NAO");
    for (int ispin = 0; ispin != mf.get_n_spins(); ispin++)
    {
        for (int isoc1 = 0; isoc1 != mf.get_n_soc(); isoc1++)
        {
            for (int isoc2 = 0; isoc2 != mf.get_n_soc(); isoc2++)
            {
                for (auto iomega = 0; iomega != tfg.get_n_grids(); iomega++)
                {
                    size_t n_IJR_myid = 0;
                    std::stringstream ss;
                    ss << input_dir << "SigcRF_ispin_" << std::setfill('0') << std::setw(2) << ispin
                       << "_s_" << std::setw(1) << isoc1 << std::setw(1) << isoc2 << "_iomega_"
                       << std::setfill('0') << std::setw(3) << iomega << "_myid_"
                       << std::setfill('0') << std::setw(5) << envs::myid_global << ".dat";
                    std::ifstream ifs_sigmac_r(ss.str(), std::ios::in | std::ios::binary);
                    if (!ifs_sigmac_r.is_open())
                    {
                        throw std::runtime_error("Cannot open file " + ss.str());
                    }
                    ifs_sigmac_r.read((char *)&n_IJR_myid, sizeof(size_t));

                    const auto omega = tfg.get_freq_nodes()[iomega];
                    for (size_t idx = 0; idx != n_IJR_myid; idx++)
                    {
                        size_t dims[5];
                        ifs_sigmac_r.read((char *)dims, 5 * sizeof(size_t));
                        const auto iR = dims[0];
                        const auto I = dims[1];
                        const auto J = dims[2];
                        const auto n_I = dims[3];
                        const auto n_J = dims[4];
                        Matz sigc(n_I, n_J, MAJOR::ROW);
                        ifs_sigmac_r.read((char *)sigc.ptr(), n_I * n_J * 2 * sizeof(double));
                        const auto R = Rlist[iR];
                        sigc_is_f_R_IJ[ispin][isoc1][isoc2][omega][R][I][J] = std::move(sigc);
                    }
                    ifs_sigmac_r.close();
                }
            }
        }
    }

    is_rspace_built_ = true;
    mpi_comm_global_h.barrier();
    LIBRPA::utils::lib_printf_root("Finished reading real-space imaginary-frequency NAO sigma_c matrices.\n");
    Profiler::stop("g0w0_read_sigc(R,iw) in NAO");
}

void G0W0::build_sigc_matrix_KS(
    const std::vector<std::vector<std::vector<ComplexMatrix>>> &wfc_target,
    const std::vector<Vector3_Order<double>> &kfrac_target)
{
    using LIBRPA::envs::blacs_ctxt_global_h;
    using LIBRPA::envs::mpi_comm_global_h;

    assert(this->is_rspace_built_);
    if (this->is_kspace_built_)
    {
        utils::lib_printf("Warning: reset Sigmac_c k-space matrices\n");
        this->reset_kspace();
    }

    const int n_aos = mf.get_n_aos();
    const int n_bands = mf.get_n_bands();

    profiler.start("g0w0_build_sigc_KS");

#ifndef LIBRPA_USE_LIBRI
    if (mpi_comm_global_h.myid == 0)
    {
        cout << "LIBRA::G0W0::build_sigc_matrix_KS is only implemented on top of LibRI" << endl;
        cout << "Please recompile LibRPA with -DUSE_LIBRI and optionally configure include path"
             << endl;
    }
    mpi_comm_global_h.barrier();
    throw std::logic_error("compilation");
#else
    // char fn[80];
    Array_Desc desc_nband_nao(blacs_ctxt_global_h);
    desc_nband_nao.init_1b1p(n_bands, n_aos, 0, 0);
    Array_Desc desc_nao_nao(blacs_ctxt_global_h);
    desc_nao_nao.init_1b1p(n_aos, n_aos, 0, 0);
    Array_Desc desc_nband_nband(blacs_ctxt_global_h);
    desc_nband_nband.init_1b1p(n_bands, n_bands, 0, 0);
    Array_Desc desc_nband_nband_fb(blacs_ctxt_global_h);
    desc_nband_nband_fb.init(n_bands, n_bands, n_bands, n_bands, 0, 0);

    // local 2D-block submatrices
    auto sigc_nao_nao = init_local_mat<complex<double>>(desc_nao_nao, MAJOR::COL);
    auto sigc_nband_nband = init_local_mat<complex<double>>(desc_nband_nband, MAJOR::COL);

    const auto set_IJ_nao_nao = LIBRPA::utils::get_necessary_IJ_from_block_2D(
        atomic_basis_wfc, atomic_basis_wfc, desc_nao_nao);
    const auto s0_s1 = get_s0_s1_for_comm_map2_first(set_IJ_nao_nao);

    // NOTE: With many MPI tasks, matrix at a particular k point
    // can have contributions from many tasks. The same happens for EXX.
    for (int ispin = 0; ispin < this->mf.get_n_spins(); ispin++)
    {
        for (int isoc1 = 0; isoc1 < this->mf.get_n_soc(); isoc1++)
        {
            for (int isoc2 = 0; isoc2 < this->mf.get_n_soc(); isoc2++)
            {
                for (const auto &freq : this->tfg.get_freq_nodes())
                {
                    const auto &sigc_is_freq =
                        this->sigc_is_f_R_IJ.at(ispin).at(isoc1).at(isoc2).at(freq);
                    // Communicate to obtain sub-matrices for necessary I-J pairs at all Rs
                    std::map<int,
                             std::map<std::pair<int, std::array<int, 3>>, Tensor<complex<double>>>>
                        sigc_I_JR_local;
                    if (this->sigc_is_f_R_IJ.count(ispin) &&
                        this->sigc_is_f_R_IJ.at(ispin).count(isoc1) &&
                        this->sigc_is_f_R_IJ.at(ispin).at(isoc1).count(isoc2))
                    {
                        const auto R = R_IJ_sigc.first;
                        for (const auto &I_J_sigc : R_IJ_sigc.second)
                        {
                            const auto I = I_J_sigc.first;
                            const auto &n_I = atomic_basis_wfc.get_atom_nb(I);
                            for (const auto &J_sigc : I_J_sigc.second)
                            {
                                const auto J = J_sigc.first;
                                const auto &n_J = atomic_basis_wfc.get_atom_nb(J);
                                const std::array<int, 3> Ra{R.x, R.y, R.z};
                                sigc_I_JR_local[I][{J, Ra}] =
                                    Tensor<complex<double>>({n_I, n_J}, J_sigc.second.sptr());
                            }
                        }
                    }
                    auto sigc_I_JR = comm_map2_first(mpi_comm_global_h.comm, sigc_I_JR_local,
                                                     s0_s1.first, s0_s1.second);
                    sigc_I_JR_local.clear();

                    // Convert each <I,<J, R>> pair to the nearest neighbour to speed up later
                    // Fourier transform while keep the accuracy in further band interpolation.
                    // Reuse the cleared-up sigc_I_JR_local object
                    if (coord_frac.size() > 0)
                    {
                        for (auto &I_sigcJR : sigc_I_JR)
                        {
                            const auto &I = I_sigcJR.first;
                            for (auto &JR_sigc : I_sigcJR.second)
                            {
                                const auto &J = JR_sigc.first.first;
                                const auto &R = JR_sigc.first.second;

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
                                                 Vector3<double>(R_IJ.x, R_IJ.y, R_IJ.z)) *
                                                latvec;
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
                                sigc_I_JR_local[I][{J, R_bvk}] = std::move(JR_sigc.second);
                            }
                        }
                    }
                    else
                    {
                        sigc_I_JR_local = std::move(sigc_I_JR);
                    }

                    // Perform Fourier transform
                    for (int ik = 0; ik < kfrac_target.size(); ik++)
                    {
                        const auto kfrac = kfrac_target[ik];

                        const std::function<complex<double>(
                            const int &, const std::pair<int, std::array<int, 3>> &)>
                            fourier = [kfrac](const int &I,
                                              const std::pair<int, std::array<int, 3>> &J_Ra)
                        {
                            const auto &Ra = J_Ra.second;
                            Vector3<double> R_IJ(Ra[0], Ra[1], Ra[2]);
                            const auto ang = (kfrac * R_IJ) * TWO_PI;
                            return complex<double>{std::cos(ang), std::sin(ang)};
                        };

                        sigc_nao_nao.zero_out();
                        collect_block_from_IJ_storage_tensor_transform(
                            sigc_nao_nao, desc_nao_nao, atomic_basis_wfc, atomic_basis_wfc, fourier,
                            sigc_I_JR_local);
                        // prepare wave function BLACS
                        const auto &wfc_isp1_k = wfc_target[ispin][isoc1][ik];
                        const auto &wfc_isp2_k = wfc_target[ispin][isoc2][ik];
                        blacs_ctxt_global_h.barrier();
                        const auto wfc1_block =
                            get_local_mat(wfc_isp1_k.c, MAJOR::ROW, desc_nband_nao, MAJOR::COL)
                                .conj();
                        const auto wfc2_block =
                            get_local_mat(wfc_isp2_k.c, MAJOR::ROW, desc_nband_nao, MAJOR::COL)
                                .conj();
                        auto temp_nband_nao = multiply_scalapack(
                            wfc1_block, desc_nband_nao, sigc_nao_nao, desc_nao_nao, desc_nband_nao);
                        ScalapackConnector::pgemm_f(
                            'N', 'C', n_bands, n_bands, n_aos, 1.0, temp_nband_nao.ptr(), 1, 1,
                            desc_nband_nao.desc, wfc2_block.ptr(), 1, 1, desc_nband_nao.desc, 0.0,
                            sigc_nband_nband.ptr(), 1, 1, desc_nband_nband.desc);
                        // collect the full matrix to master
                        // TODO: would need a different strategy for large system
                        auto sigc_nband_nband_fb =
                            init_local_mat<complex<double>>(desc_nband_nband_fb, MAJOR::COL);
                        ScalapackConnector::pgemr2d_f(
                            n_bands, n_bands, sigc_nband_nband.ptr(), 1, 1, desc_nband_nband.desc,
                            sigc_nband_nband_fb.ptr(), 1, 1, desc_nband_nband_fb.desc,
                            desc_nband_nband_fb.ictxt());
                        // NOTE: only the matrices at master process is meaningful
                        if (sigc_is_ik_f_KS.count(ispin) == 0 ||
                            sigc_is_ik_f_KS.at(ispin).count(ik) == 0 ||
                            sigc_is_ik_f_KS.at(ispin).at(ik).count(freq) == 0)
                        {
                            sigc_is_ik_f_KS[ispin][ik][freq] =
                                init_local_mat<complex<double>>(desc_nband_nband_fb, MAJOR::COL);
                        }
                        sigc_is_ik_f_KS[ispin][ik][freq] += sigc_nband_nband_fb;
                    }
                }
            }
        }
        //负频
        for (const auto& freq: this->tfg.get_freq_nodes())
        {
            const auto& sigc_is_freq = this->sigc_is_f_R_IJ.at(ispin).at(-freq);
            // Communicate to obtain sub-matrices for necessary I-J pairs at all Rs
            std::map<int, std::map<std::pair<int, std::array<int, 3>>, Tensor<complex<double>>>>
                sigc_I_JR_local;
            for (const auto &R_IJ_sigc: sigc_is_freq)
            {
                const auto R = R_IJ_sigc.first;
                for (const auto &I_J_sigc: R_IJ_sigc.second)
                {
                    const auto I = I_J_sigc.first;
                    const auto &n_I = atomic_basis_wfc.get_atom_nb(I);
                    for (const auto &J_sigc: I_J_sigc.second)
                    {
                        const auto J = J_sigc.first;
                        const auto &n_J = atomic_basis_wfc.get_atom_nb(J);
                        const std::array<int, 3> Ra{R.x, R.y, R.z};
                        sigc_I_JR_local[I][{J, Ra}] = Tensor<complex<double>>({n_I, n_J}, J_sigc.second.sptr());
                    }
                }
            }
            auto sigc_I_JR = comm_map2_first(mpi_comm_global_h.comm, sigc_I_JR_local, s0_s1.first, s0_s1.second);
            sigc_I_JR_local.clear();
            
            // Convert each <I,<J, R>> pair to the nearest neighbour to speed up later Fourier transform
            // while keep the accuracy in further band interpolation.
            // Reuse the cleared-up sigc_I_JR_local object
            if (coord_frac.size() > 0)
            {
                for (auto &I_sigcJR: sigc_I_JR)
                {
                    const auto &I = I_sigcJR.first;
                    for (auto &JR_sigc: I_sigcJR.second)
                    {
                        const auto &J = JR_sigc.first.first;
                        const auto &R = JR_sigc.first.second;

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
                        sigc_I_JR_local[I][{J, R_bvk}] = std::move(JR_sigc.second);
                    }
                }
            }
            else
            {
                sigc_I_JR_local = std::move(sigc_I_JR);
            }

            // Perform Fourier transform
            for (int ik = 0; ik < kfrac_target.size(); ik++)
            {
                const auto kfrac = kfrac_target[ik];

                const std::function<complex<double>(const int &, const std::pair<int, std::array<int, 3>> &)>
                    fourier = [kfrac](const int &I, const std::pair<int, std::array<int, 3>> &J_Ra)
                    {
                        const auto &Ra = J_Ra.second;
                        Vector3<double> R_IJ(Ra[0], Ra[1], Ra[2]);
                        const auto ang = (kfrac * R_IJ) * TWO_PI;
                        return complex<double>{std::cos(ang), std::sin(ang)};
                    };

                sigc_nao_nao.zero_out();
                collect_block_from_IJ_storage_tensor_transform(sigc_nao_nao, desc_nao_nao, 
                        atomic_basis_wfc, atomic_basis_wfc,
                        fourier, sigc_I_JR_local);
                // prepare wave function BLACS
                const auto &wfc_isp_k = wfc_target[ispin][ik];
                blacs_ctxt_global_h.barrier();
                const auto wfc_block = get_local_mat(wfc_isp_k.c, MAJOR::ROW, desc_nband_nao, MAJOR::COL).conj();
                auto temp_nband_nao = multiply_scalapack(wfc_block, desc_nband_nao, sigc_nao_nao, desc_nao_nao, desc_nband_nao);
                ScalapackConnector::pgemm_f('N', 'C', n_bands, n_bands, n_aos, 1.0,
                                            temp_nband_nao.ptr(), 1, 1, desc_nband_nao.desc,
                                            wfc_block.ptr(), 1, 1, desc_nband_nao.desc, 0.0,
                                            sigc_nband_nband.ptr(), 1, 1, desc_nband_nband.desc);
                // collect the full matrix to master
                // TODO: would need a different strategy for large system
                auto sigc_nband_nband_fb = init_local_mat<complex<double>>(desc_nband_nband_fb, MAJOR::COL);
                ScalapackConnector::pgemr2d_f(n_bands, n_bands,
                                              sigc_nband_nband.ptr(), 1, 1, desc_nband_nband.desc,
                                              sigc_nband_nband_fb.ptr(), 1, 1, desc_nband_nband_fb.desc,
                                              desc_nband_nband_fb.ictxt());
                // NOTE: only the matrices at master process is meaningful
                sigc_is_ik_f_KS[ispin][ik][-freq] = sigc_nband_nband_fb;
            }
        }
    }
#endif
    this->is_kspace_built_ = true;
    profiler.stop("g0w0_build_sigc_KS");
}

void G0W0::build_sigc_matrix_KS_kgrid()
{
    LIBRPA::envs::mpi_comm_global_h.barrier();
    if (LIBRPA::envs::mpi_comm_global_h.myid == 0)
    {
        LIBRPA::utils::lib_printf(
            "build_sigc_matrix_KS_kgrid: constructing self-energy matrix for SCF k-grid\n");
    }
    this->build_sigc_matrix_KS(this->mf.get_eigenvectors(), this->kfrac_list);
}

void G0W0::build_sigc_matrix_KS_kgrid0()
{
    this->build_sigc_matrix_KS(this->mf.get_eigenvectors0(), this->kfrac_list);
}

void G0W0::build_sigc_matrix_KS_band(
    const std::vector<std::vector<std::vector<ComplexMatrix>>> &wfc,
    const std::vector<Vector3_Order<double>> &kfrac_band)
{
    LIBRPA::envs::mpi_comm_global_h.barrier();
    if (LIBRPA::envs::mpi_comm_global_h.myid == 0)
    {
        LIBRPA::utils::lib_printf(
            "build_sigc_matrix_KS_kgrid: constructing self-energy matrix for band k-path\n");
    }
    this->build_sigc_matrix_KS(wfc, kfrac_band);
}

template void G0W0::build_spacetime<double>(
    const Cs_LRI &,
    map<double,
        atom_mapping<std::map<Vector3_Order<double>, matrix_m<complex<double>>>>::pair_t_old> &,
    const vector<Vector3_Order<int>> &, const vector<Vector3_Order<double>> &,
    std::map<Vector3_Order<double>, ComplexMatrix> &);
template void G0W0::build_spacetime<std::complex<double>>(
    const Cs_LRI &,
    map<double,
        atom_mapping<std::map<Vector3_Order<double>, matrix_m<complex<double>>>>::pair_t_old> &,
    const vector<Vector3_Order<int>> &, const vector<Vector3_Order<double>> &,
    std::map<Vector3_Order<double>, ComplexMatrix> &);

}  // namespace LIBRPA
