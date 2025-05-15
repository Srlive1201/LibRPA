#include "task_gw_band.h"

#include <fstream>
#include <sstream>

#include "analycont.h"
#include "chi0.h"
#include "constants.h"
#include "coulmat.h"
#include "driver_params.h"
#include "driver_utils.h"
#include "envs_io.h"
#include "envs_mpi.h"
#include "epsilon.h"
#include "exx.h"
#include "gw.h"
#include "meanfield.h"
#include "params.h"
#include "pbc.h"
#include "profiler.h"
#include "qpe_solver.h"
#include "read_data.h"
#include "ri.h"
#include "utils_timefreq.h"
#include "write_aims.h"

void task_g0w0_band()
{
    using LIBRPA::envs::mpi_comm_global_h;
    using LIBRPA::envs::ofs_myid;
    using LIBRPA::utils::lib_printf;

    Profiler::start("g0w0_band", "G0W0 quasi-particle band structure calculation");

    Vector3_Order<int> period{kv_nmp[0], kv_nmp[1], kv_nmp[2]};
    auto Rlist = construct_R_grid(period);

    vector<Vector3_Order<double>> qlist;
    for (auto q_weight : irk_weight)
    {
        qlist.push_back(q_weight.first);
    }

    // Prepare time-frequency grids
    auto tfg =
        LIBRPA::utils::generate_timefreq_grids(Params::nfreq, Params::tfgrids_type, meanfield);

    Chi0 chi0(meanfield, klist, tfg);
    chi0.gf_R_threshold = Params::gf_R_threshold;

    Profiler::start("chi0_build", "Build response function chi0");
    if (Params::use_shrink_abfs)
    {
        // replace atom_mu by atom_mu_l to construct chi0 due to LRI error
        atom_mu = atom_mu_l;
        LIBRPA::atomic_basis_abf.set(atom_mu);
        atom_mu_part_range.resize(atom_mu.size());
        atom_mu_part_range[0] = 0;
        for (int I = 1; I != atom_mu.size(); I++)
            atom_mu_part_range[I] = atom_mu.at(I - 1) + atom_mu_part_range[I - 1];

        N_all_mu = atom_mu_part_range[natom - 1] + atom_mu[natom - 1];
    }
    chi0.build(Cs_data, Rlist, period, local_atpair, qlist);
    Profiler::stop("chi0_build");

    // for (auto &Ip : chi0.get_chi0_q().at(tfg.get_freq_nodes()[10]).at(qlist[0]))
    // {
    //     auto I = Ip.first;
    //     for (auto &Jm : Ip.second)
    //     {
    //         auto J = Jm.first;
    //         print_complex_matrix_mm(Jm.second,
    //                                 "chi0_" + std::to_string(I) + "_" + std::to_string(J), 0.);
    //     }
    // }

    std::flush(ofs_myid);
    mpi_comm_global_h.barrier();

    std::map<Vector3_Order<double>, ComplexMatrix> sinvS;
    if (Params::use_shrink_abfs)
    {
        Profiler::start("read_shrink_sinvS_fold", "Load shrink transformation");
        // change atom_mu: number of {Mu,mu} in the later calculations
        read_shrink_sinvS(driver_params.input_dir, "shrink_sinvS_", sinvS);
        Profiler::stop("read_shrink_sinvS_fold");
        Profiler::start("shrink_chi0_abfs", "Do shrink transformation");
        chi0.shrink_abfs_chi0(sinvS, qlist, atom_mu_l, atom_mu);
        Profiler::stop("shrink_chi0_abfs");
    }

    if (Params::debug)
    {  // debug, check chi0
        char fn[80];
        for (const auto &chi0q : chi0.get_chi0_q())
        {
            const int ifreq = chi0.tfg.get_freq_index(chi0q.first);
            for (const auto &q_IJchi0 : chi0q.second)
            {
                const int iq = std::distance(klist.begin(),
                                             std::find(klist.begin(), klist.end(), q_IJchi0.first));
                for (const auto &I_Jchi0 : q_IJchi0.second)
                {
                    const auto &I = I_Jchi0.first;
                    for (const auto &J_chi0 : I_Jchi0.second)
                    {
                        const auto &J = J_chi0.first;
                        sprintf(fn, "chi0fq_ifreq_%d_iq_%d_I_%zu_J_%zu_id_%d.mtx", ifreq, iq, I, J,
                                mpi_comm_global_h.myid);
                        print_complex_matrix_mm(J_chi0.second, Params::output_dir + "/" + fn,
                                                1e-15);
                    }
                }
            }
        }
    }

    Profiler::start("read_vq_cut", "Load truncated Coulomb");
    if (LIBRPA::parallel_routing == LIBRPA::ParallelRouting::R_TAU)
    {
        read_Vq_full(driver_params.input_dir, "coulomb_cut_", true);
    }
    else
    {
        // NOTE: local_atpair already set in the main.cpp.
        //       It can consists of distributed atom pairs of only upper half.
        //       Setup of local_atpair may be better to extracted as some util function,
        //       instead of in the main driver.
        read_Vq_row(driver_params.input_dir, "coulomb_cut_", Params::vq_threshold, local_atpair,
                    true);
    }
    Profiler::cease("read_vq_cut");

    std::vector<double> epsmac_LF_imagfreq_re;
    if (Params::replace_w_head)
    {
        std::vector<double> omegas_dielect;
        std::vector<double> dielect_func;
        if (Params::option_dielect_func != 3 && Params::option_dielect_func != 4)
            read_dielec_func(driver_params.input_dir + "dielecfunc_out", omegas_dielect,
                             dielect_func);

        epsmac_LF_imagfreq_re = interpolate_dielec_func(Params::option_dielect_func, omegas_dielect,
                                                        dielect_func, chi0.tfg.get_freq_nodes());
    }

    Profiler::start("g0w0_exx", "Build exchange self-energy");
    auto exx = LIBRPA::Exx(meanfield, kfrac_list, period);
    {
        Profiler::start("ft_vq_cut", "Fourier transform truncated Coulomb");
        const auto VR = FT_Vq(Vq_cut, meanfield.get_n_kpoints(), Rlist, true);
        Profiler::stop("ft_vq_cut");

        Profiler::start("g0w0_exx_real_work");
        if (Params::use_shrink_abfs)
        {
            if (Params::use_soc)
                exx.build<std::complex<double>>(Cs_shrinked_data, Rlist, VR);
            else
                exx.build<double>(Cs_shrinked_data, Rlist, VR);
        }
        else
        {
            if (Params::use_soc)
                exx.build<std::complex<double>>(Cs_data, Rlist, VR);
            else
                exx.build<double>(Cs_data, Rlist, VR);
        }
        Profiler::stop("g0w0_exx_real_work");
    }
    Profiler::stop("g0w0_exx");
    std::flush(ofs_myid);

    Profiler::start("g0w0_wc", "Build screened interaction");
    vector<std::complex<double>> epsmac_LF_imagfreq(epsmac_LF_imagfreq_re.cbegin(),
                                                    epsmac_LF_imagfreq_re.cend());
    map<double,
        atom_mapping<std::map<Vector3_Order<double>, matrix_m<complex<double>>>>::pair_t_old>
        Wc_freq_q;
    if (Params::use_scalapack_gw_wc)
    {
        Wc_freq_q = compute_Wc_freq_q_blacs(chi0, Vq, Vq_cut, epsmac_LF_imagfreq);
    }
    else
    {
        Wc_freq_q = compute_Wc_freq_q(chi0, Vq, Vq_cut, epsmac_LF_imagfreq);
    }
    Profiler::stop("g0w0_wc");
    if (Params::debug)
    {  // debug, check Wc
        char fn[80];
        for (const auto &Wc : Wc_freq_q)
        {
            const int ifreq = chi0.tfg.get_freq_index(Wc.first);
            for (const auto &I_JqWc : Wc.second)
            {
                const auto &I = I_JqWc.first;
                for (const auto &J_qWc : I_JqWc.second)
                {
                    const auto &J = J_qWc.first;
                    for (const auto &q_Wc : J_qWc.second)
                    {
                        const int iq = std::distance(
                            klist.begin(), std::find(klist.begin(), klist.end(), q_Wc.first));
                        sprintf(fn, "Wcfq_ifreq_%d_iq_%d_I_%zu_J_%zu_id_%d.mtx", ifreq, iq, I, J,
                                mpi_comm_global_h.myid);
                        print_matrix_mm_file(q_Wc.second, Params::output_dir + "/" + fn, 1e-15);
                    }
                }
            }
        }
    }

    if (Params::use_shrink_abfs)
    {
        auto atom_mu_s = atom_mu;
        atom_mu = atom_mu_l;
        LIBRPA::atomic_basis_abf.set(atom_mu);
        atom_mu_part_range.resize(atom_mu.size());
        atom_mu_part_range[0] = 0;
        for (int I = 1; I != atom_mu.size(); I++)
            atom_mu_part_range[I] = atom_mu.at(I - 1) + atom_mu_part_range[I - 1];

        N_all_mu = atom_mu_part_range[natom - 1] + atom_mu[natom - 1];
        Profiler::start("unfold_Wc_abfs", "Do shrink transformation");
        chi0.unfold_abfs_Wc(sinvS, Wc_freq_q, qlist, atom_mu, atom_mu_s);
        Profiler::stop("unfold_Wc_abfs");
        sinvS.clear();
    }

    LIBRPA::G0W0 s_g0w0(meanfield, kfrac_list, chi0.tfg, period);
    Profiler::start("g0w0_sigc_IJ", "Build real-space correlation self-energy");

    if (Params::use_soc)
        s_g0w0.build_spacetime<std::complex<double>>(Cs_data, Wc_freq_q, Rlist);
    else
        s_g0w0.build_spacetime<double>(Cs_data, Wc_freq_q, Rlist);

    Profiler::stop("g0w0_sigc_IJ");
    std::flush(ofs_myid);

    /*
     * Compute the QP energies on k-grid
     */
    Profiler::start("read_vxc", "Load DFT xc potential");
    std::vector<matrix> vxc;
    int flag_read_vxc = read_vxc(driver_params.input_dir + "vxc_out", vxc);
    Profiler::stop("read_vxc");

    if (flag_read_vxc != 0)
    {
        if (mpi_comm_global_h.myid == 0)
        {
            lib_printf("Error in reading Vxc on kgrid, task failed!\n");
        }
        Profiler::stop("g0w0_band");
        return;
    }

    Profiler::start("g0w0_exx_ks_kgrid");
    exx.build_KS_kgrid();
    Profiler::stop("g0w0_exx_ks_kgrid");
    Profiler::start("g0w0_sigc_ks_kgrid");
    s_g0w0.build_sigc_matrix_KS_kgrid();
    Profiler::stop("g0w0_sigc_ks_kgrid");

    // imaginary freqencies for analytic continuation
    std::vector<cplxdb> imagfreqs;
    for (const auto &freq : chi0.tfg.get_freq_nodes())
    {
        imagfreqs.push_back(cplxdb{0.0, freq});
    }

    Profiler::start("g0w0_solve_qpe_kgrid", "Solve quasi-particle equation");
    if (mpi_comm_global_h.is_root())
    {
        std::cout << "Solving quasi-particle equation for states at k-points of regular grid\n";
    }

    // TODO: parallelize analytic continuation and QPE solver among tasks
    if (mpi_comm_global_h.is_root())
    {
        map<int, map<int, map<int, double>>> e_qp_all;
        map<int, map<int, map<int, cplxdb>>> sigc_all;
        const auto efermi = meanfield.get_efermi();
        for (int i_spin = 0; i_spin < meanfield.get_n_spins(); i_spin++)
        {
            for (int i_kpoint = 0; i_kpoint < meanfield.get_n_kpoints(); i_kpoint++)
            {
                const auto &sigc_sk = s_g0w0.sigc_is_ik_f_KS[i_spin][i_kpoint];
                for (int i_state = 0; i_state < meanfield.get_n_bands(); i_state++)
                {
                    const auto &eks_state = meanfield.get_eigenvals()[i_spin](i_kpoint, i_state);
                    const auto &exx_state = exx.Eexx[i_spin][i_kpoint][i_state];
                    const auto &vxc_state = vxc[i_spin](i_kpoint, i_state);
                    std::vector<cplxdb> sigc_state;
                    for (const auto &freq : chi0.tfg.get_freq_nodes())
                    {
                        sigc_state.push_back(sigc_sk.at(freq)(i_state, i_state));
                    }
                    LIBRPA::AnalyContPade pade(Params::n_params_anacon, imagfreqs, sigc_state);
                    double e_qp;
                    cplxdb sigc;
                    int flag_qpe_solver = LIBRPA::qpe_solver_pade_self_consistent(
                        pade, eks_state, efermi, vxc_state, exx_state, e_qp, sigc);
                    if (flag_qpe_solver == 0)
                    {
                        e_qp_all[i_spin][i_kpoint][i_state] = e_qp;
                        sigc_all[i_spin][i_kpoint][i_state] = sigc;
                    }
                    else
                    {
                        printf("Warning! QPE solver failed for spin %d, kpoint %d, state %d\n",
                               i_spin + 1, i_kpoint + 1, i_state + 1);
                        e_qp_all[i_spin][i_kpoint][i_state] =
                            std::numeric_limits<double>::quiet_NaN();
                        sigc_all[i_spin][i_kpoint][i_state] =
                            std::numeric_limits<cplxdb>::quiet_NaN();
                    }
                }
            }
        }

        // display results
        const std::string banner(124, '-');
        printf("Printing quasi-particle energy [unit: eV]\n\n");
        for (int i_spin = 0; i_spin < meanfield.get_n_spins(); i_spin++)
        {
            for (int i_kpoint = 0; i_kpoint < meanfield.get_n_kpoints(); i_kpoint++)
            {
                const auto &k = kfrac_list[i_kpoint];
                printf("spin %2d, k-point %4d: (%.5f, %.5f, %.5f) \n", i_spin + 1, i_kpoint + 1,
                       k.x, k.y, k.z);
                printf("%124s\n", banner.c_str());
                printf("%5s %16s %16s %16s %16s %16s %16s %16s\n", "State", "occ", "e_mf", "v_xc",
                       "v_exx", "ReSigc", "ImSigc", "e_qp");
                printf("%124s\n", banner.c_str());
                for (int i_state = 0; i_state < meanfield.get_n_bands(); i_state++)
                {
                    const auto &occ_state = meanfield.get_weight()[i_spin](i_kpoint, i_state) *
                                            meanfield.get_n_kpoints();
                    const auto &eks_state =
                        meanfield.get_eigenvals()[i_spin](i_kpoint, i_state) * HA2EV;
                    const auto &exx_state = exx.Eexx[i_spin][i_kpoint][i_state] * HA2EV;
                    const auto &vxc_state = vxc[i_spin](i_kpoint, i_state) * HA2EV;
                    const auto &resigc = sigc_all[i_spin][i_kpoint][i_state].real() * HA2EV;
                    const auto &imsigc = sigc_all[i_spin][i_kpoint][i_state].imag() * HA2EV;
                    const auto &eqp = e_qp_all[i_spin][i_kpoint][i_state] * HA2EV;
                    printf("%5d %16.5f %16.5f %16.5f %16.5f %16.5f %16.5f %16.5f\n", i_state + 1,
                           occ_state, eks_state, vxc_state, exx_state, resigc, imsigc, eqp);
                }
                printf("\n");
            }
        }
    }
    Profiler::stop("g0w0_solve_qpe_kgrid");

    /*
     * Compute the QP energies on band k-paths
     */
    // Reset k-space EXX and Sigmac matrices to avoid warning from internal reset
    exx.reset_kspace();
    s_g0w0.reset_kspace();

    /* Below we handle the band k-points data
     * First load the information of k-points along the k-path */
    Profiler::start("g0w0_band_load_band_mf", "Read eigen solutions at band kpoints");
    int n_basis_band, n_states_band, n_spin_band;
    std::vector<Vector3_Order<double>> kfrac_band = read_band_kpath_info(
        driver_params.input_dir + "band_kpath_info", n_basis_band, n_states_band, n_spin_band);
    if (mpi_comm_global_h.is_root())
    {
        std::cout << "Band k-points to compute:\n";
        for (int ik = 0; ik < kfrac_band.size(); ik++)
        {
            const auto &k = kfrac_band[ik];
            lib_printf("%5d %12.7f %12.7f %12.7f\n", ik + 1, k.x, k.y, k.z);
        }
    }
    mpi_comm_global_h.barrier();

    auto meanfield_band = read_meanfield_band(driver_params.input_dir, n_basis_band, n_states_band,
                                              n_spin_band, kfrac_band.size());

    /* Set the same Fermi energy as in SCF */
    meanfield_band.get_efermi() = meanfield.get_efermi();
    Profiler::stop("g0w0_band_load_band_mf");

    Profiler::start("g0w0_sigx_rotate_KS");
    exx.build_KS_band(meanfield_band.get_eigenvectors(), kfrac_band);
    Profiler::stop("g0w0_sigx_rotate_KS");
    std::flush(ofs_myid);

    Profiler::start("g0w0_sigc_rotate_KS");
    s_g0w0.build_sigc_matrix_KS_band(meanfield_band.get_eigenvectors(), kfrac_band);
    Profiler::stop("g0w0_sigc_rotate_KS");
    std::flush(ofs_myid);

    Profiler::start("read_vxc", "Load DFT xc potential");
    auto vxc_band =
        read_vxc_band(driver_params.input_dir, n_states_band, n_spin_band, kfrac_band.size());
    Profiler::stop("read_vxc");
    std::flush(ofs_myid);

    Profiler::start("g0w0_solve_qpe", "Solve quasi-particle equation");
    if (mpi_comm_global_h.is_root())
    {
        std::cout << "Solving quasi-particle equation\n";
    }

    // TODO: parallelize analytic continuation and QPE solver among tasks
    if (mpi_comm_global_h.is_root())
    {
        const auto &mf = meanfield_band;
        map<int, map<int, map<int, double>>> e_qp_all;
        map<int, map<int, map<int, cplxdb>>> sigc_all;
        const auto efermi = mf.get_efermi();
        for (int i_spin = 0; i_spin < mf.get_n_spins(); i_spin++)
        {
            for (int i_kpoint = 0; i_kpoint < mf.get_n_kpoints(); i_kpoint++)
            {
                const auto &sigc_sk = s_g0w0.sigc_is_ik_f_KS[i_spin][i_kpoint];
                for (int i_state = 0; i_state < mf.get_n_bands(); i_state++)
                {
                    const auto &eks_state = mf.get_eigenvals()[i_spin](i_kpoint, i_state);
                    const auto &exx_state = exx.Eexx[i_spin][i_kpoint][i_state];
                    const auto &vxc_state = vxc_band[i_spin](i_kpoint, i_state);
                    std::vector<cplxdb> sigc_state;
                    for (const auto &freq : chi0.tfg.get_freq_nodes())
                    {
                        sigc_state.push_back(sigc_sk.at(freq)(i_state, i_state));
                    }
                    LIBRPA::AnalyContPade pade(Params::n_params_anacon, imagfreqs, sigc_state);
                    double e_qp;
                    cplxdb sigc;
                    int flag_qpe_solver = LIBRPA::qpe_solver_pade_self_consistent(
                        pade, eks_state, efermi, vxc_state, exx_state, e_qp, sigc);
                    if (flag_qpe_solver == 0)
                    {
                        e_qp_all[i_spin][i_kpoint][i_state] = e_qp;
                        sigc_all[i_spin][i_kpoint][i_state] = sigc;
                    }
                    else
                    {
                        printf("Warning! QPE solver failed for spin %d, kpoint %d, state %d\n",
                               i_spin + 1, i_kpoint + 1, i_state + 1);
                        e_qp_all[i_spin][i_kpoint][i_state] =
                            std::numeric_limits<double>::quiet_NaN();
                        sigc_all[i_spin][i_kpoint][i_state] =
                            std::numeric_limits<cplxdb>::quiet_NaN();
                    }
                }
            }
        }
        // output bandgap
        double bandgap = 0.0;
        double valence = -1.e10;
        double conduct = 1.e10;
        int ik_val = 0;
        int ik_cond = 0;
        int nocc = 0;
        auto &wg = meanfield.get_weight()[0];
        for (int i = 0; i != wg.size; i++)
        {
            if (wg.c[i] == 0.)
            {
                nocc = i;
                break;
            }
        }
        lib_printf("Bands of occupation: %4d \n", nocc);

        // display results
        for (int i_spin = 0; i_spin < mf.get_n_spins(); i_spin++)
        {
            std::ofstream ofs_ks;
            std::ofstream ofs_hf;
            std::ofstream ofs_gw;
            std::stringstream fn;

            fn << "GW_band_spin_" << i_spin + 1 << ".dat";
            ofs_gw.open(fn.str());

            fn.str("");
            fn.clear();
            fn << "EXX_band_spin_" << i_spin + 1 << ".dat";
            ofs_hf.open(fn.str());

            fn.str("");
            fn.clear();
            fn << "KS_band_spin_" << i_spin + 1 << ".dat";
            ofs_ks.open(fn.str());

            ofs_gw << std::fixed;
            ofs_hf << std::fixed;
            ofs_ks << std::fixed;

            for (int i_kpoint = 0; i_kpoint < mf.get_n_kpoints(); i_kpoint++)
            {
                const auto &k = kfrac_band[i_kpoint];
                ofs_ks << std::setw(5) << i_kpoint + 1 << std::setw(15) << std::setprecision(7)
                       << k.x << std::setw(15) << std::setprecision(7) << k.y << std::setw(15)
                       << std::setprecision(7) << k.z;
                ofs_gw << std::setw(5) << i_kpoint + 1 << std::setw(15) << std::setprecision(7)
                       << k.x << std::setw(15) << std::setprecision(7) << k.y << std::setw(15)
                       << std::setprecision(7) << k.z;
                ofs_hf << std::setw(5) << i_kpoint + 1 << std::setw(15) << std::setprecision(7)
                       << k.x << std::setw(15) << std::setprecision(7) << k.y << std::setw(15)
                       << std::setprecision(7) << k.z;
                for (int i_state = 0; i_state < meanfield.get_n_bands(); i_state++)
                {
                    const auto &occ_state = mf.get_weight()[i_spin](i_kpoint, i_state);
                    const auto &eks_state = mf.get_eigenvals()[i_spin](i_kpoint, i_state) * HA2EV;
                    const auto &exx_state = exx.Eexx[i_spin][i_kpoint][i_state] * HA2EV;
                    const auto &vxc_state = vxc_band[i_spin](i_kpoint, i_state) * HA2EV;
                    // const auto &resigc = sigc_all[i_spin][i_kpoint][i_state].real() * HA2EV;
                    // const auto &imsigc = sigc_all[i_spin][i_kpoint][i_state].imag() * HA2EV;
                    const auto &eqp = e_qp_all[i_spin][i_kpoint][i_state] * HA2EV;
                    ofs_ks << std::setw(15) << std::setprecision(5) << occ_state << std::setw(15)
                           << std::setprecision(5) << eks_state;
                    ofs_gw << std::setw(15) << std::setprecision(5) << occ_state << std::setw(15)
                           << std::setprecision(5) << eqp;
                    ofs_hf << std::setw(15) << std::setprecision(5) << occ_state << std::setw(15)
                           << std::setprecision(5) << eks_state - vxc_state + exx_state;

                    // output bandgap
                    if (i_state == nocc - 1 && eqp > valence)  // HOMO
                    {
                        valence = eqp;
                        ik_val = i_kpoint;
                    }
                    else if (i_state == nocc && eqp < conduct)  // LUMO
                    {
                        conduct = eqp;
                        ik_cond = i_kpoint;
                    }
                }
                ofs_gw << "\n";
                ofs_hf << "\n";
                ofs_ks << "\n";
            }
        }
        bandgap = conduct - valence;
        const auto &k_val = kfrac_list[ik_val];
        printf("VBM: k-point %4d: (%.5f, %.5f, %.5f) \n", ik_val + 1, k_val.x, k_val.y, k_val.z);
        const auto &k_cond = kfrac_list[ik_cond];
        printf("CBM: k-point %4d: (%.5f, %.5f, %.5f) \n", ik_cond + 1, k_cond.x, k_cond.y,
               k_cond.z);
        lib_printf("Bandgap(eV): %12.7f \n", bandgap);
    }
    Profiler::stop("g0w0_solve_qpe");

    Profiler::stop("g0w0_band");
}
