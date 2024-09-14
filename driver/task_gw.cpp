#include "task_gw.h"

#include "meanfield.h"
#include "params.h"
#include "pbc.h"
#include "chi0.h"
#include "gw.h"
#include "analycont.h"
#include "qpe_solver.h"
#include "epsilon.h"
#include "exx.h"
#include "constants.h"
#include "coulmat.h"
#include "profiler.h"
#include "ri.h"

#include "envs_mpi.h"
#include "envs_io.h"
#include "utils_timefreq.h"

#include "read_data.h"
#include "write_aims.h"
#include "driver_utils.h"

void task_g0w0()
{
    using LIBRPA::envs::mpi_comm_global_h;
    using LIBRPA::envs::ofs_myid;
    using LIBRPA::utils::lib_printf;

    Profiler::start("g0w0", "G0W0 quasi-particle calculation");

    Vector3_Order<int> period {kv_nmp[0], kv_nmp[1], kv_nmp[2]};
    auto Rlist = construct_R_grid(period);

    vector<Vector3_Order<double>> qlist;
    for ( auto q_weight: irk_weight)
    {
        qlist.push_back(q_weight.first);
    }

    // Prepare time-frequency grids
    auto tfg = LIBRPA::utils::generate_timefreq_grids(Params::nfreq, Params::tfgrids_type, meanfield);

    Chi0 chi0(meanfield, klist, tfg);
    chi0.gf_R_threshold = Params::gf_R_threshold;

    Profiler::start("chi0_build", "Build response function chi0");
    chi0.build(Cs_data, Rlist, period, local_atpair, qlist);
    Profiler::stop("chi0_build");

    std::flush(ofs_myid);
    mpi_comm_global_h.barrier();

    if (Params::debug)
    { // debug, check chi0
        char fn[80];
        for (const auto &chi0q: chi0.get_chi0_q())
        {
            const int ifreq = chi0.tfg.get_freq_index(chi0q.first);
            for (const auto &q_IJchi0: chi0q.second)
            {
                const int iq = std::distance(klist.begin(), std::find(klist.begin(), klist.end(), q_IJchi0.first));
                for (const auto &I_Jchi0: q_IJchi0.second)
                {
                    const auto &I = I_Jchi0.first;
                    for (const auto &J_chi0: I_Jchi0.second)
                    {
                        const auto &J = J_chi0.first;
                        sprintf(fn, "chi0fq_ifreq_%d_iq_%d_I_%zu_J_%zu_id_%d.mtx", ifreq, iq, I, J, mpi_comm_global_h.myid);
                        print_complex_matrix_mm(J_chi0.second, Params::output_dir + "/" + fn, 1e-15);
                    }
                }
            }
        }
    }

    Profiler::start("read_vq_cut", "Load truncated Coulomb");
    read_Vq_full("./", "coulomb_cut_", true);
    Profiler::stop("read_vq_cut");

    Profiler::start("read_vxc", "Load DFT xc potential");
    std::vector<matrix> vxc;
    int flag_read_vxc = read_vxc("./vxc_out", vxc);
    if (mpi_comm_global_h.myid == 0)
    {
        if (flag_read_vxc == 0)
        {
            cout << "* Success: Read DFT xc potential, will solve quasi-particle equation\n";
        }
        else
        {
            cout << "*   Error: Read DFT xc potential, switch off solving quasi-particle equation\n";
        }
    }
    Profiler::stop("read_vxc");
    std::flush(ofs_myid);

    std::vector<double> epsmac_LF_imagfreq_re;

    if (Params::replace_w_head)
    {
        std::vector<double> omegas_dielect;
        std::vector<double> dielect_func;
        read_dielec_func("dielecfunc_out", omegas_dielect, dielect_func);

        epsmac_LF_imagfreq_re = interpolate_dielec_func(
                Params::option_dielect_func, omegas_dielect, dielect_func,
                chi0.tfg.get_freq_nodes());

        if (Params::debug)
        {
            if (mpi_comm_global_h.is_root())
            {
                lib_printf("Dielectric function parsed:\n");
                for (int i = 0; i < chi0.tfg.get_freq_nodes().size(); i++)
                    lib_printf("%d %f %f\n", i+1, chi0.tfg.get_freq_nodes()[i], epsmac_LF_imagfreq_re[i]);
            }
            mpi_comm_global_h.barrier();
        }
    }

    Profiler::start("g0w0_exx", "Build exchange self-energy");
    auto exx = LIBRPA::Exx(meanfield, kfrac_list);
    {
        Profiler::start("ft_vq_cut", "Fourier transform truncated Coulomb");
        const auto VR = FT_Vq(Vq_cut, Rlist, true);
        Profiler::stop("ft_vq_cut");

        Profiler::start("g0w0_exx_real_work");
        exx.build(Cs_data, Rlist, period, VR);
        exx.build_KS_kgrid();
        Profiler::stop("g0w0_exx_real_work");
    }
    Profiler::stop("g0w0_exx");
    std::flush(ofs_myid);

    // Since EXX contribution will be printed along with sigma_c when Vxc is read
    // we print it here when Vxc read failed
    if (flag_read_vxc != 0 && mpi_comm_global_h.is_root())
    {
        for (int isp = 0; isp != meanfield.get_n_spins(); isp++)
        {
            lib_printf("Spin channel %1d\n", isp+1);
            for (int ik = 0; ik != meanfield.get_n_kpoints(); ik++)
            {
                cout << "k-point " << ik + 1 << ": " << kvec_c[ik] << endl;
                lib_printf("%-4s  %-10s  %-10s\n", "Band", "e_exx (Ha)", "e_exx (eV)");
                for (int ib = 0; ib != meanfield.get_n_bands(); ib++)
                    lib_printf("%4d  %10.5f  %10.5f\n", ib+1, exx.Eexx[isp][ik][ib], HA2EV * exx.Eexx[isp][ik][ib]);
                lib_printf("\n");
            }
        }
    }
    mpi_comm_global_h.barrier();

    Profiler::start("g0w0_wc", "Build screened interaction");
    vector<std::complex<double>> epsmac_LF_imagfreq(epsmac_LF_imagfreq_re.cbegin(), epsmac_LF_imagfreq_re.cend());
    map<double, atom_mapping<std::map<Vector3_Order<double>, matrix_m<complex<double>>>>::pair_t_old> Wc_freq_q;
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
    { // debug, check Wc
        char fn[80];
        for (const auto &Wc: Wc_freq_q)
        {
            const int ifreq = chi0.tfg.get_freq_index(Wc.first);
            for (const auto &I_JqWc: Wc.second)
            {
                const auto &I = I_JqWc.first;
                for (const auto &J_qWc: I_JqWc.second)
                {
                    const auto &J = J_qWc.first;
                    for (const auto &q_Wc: J_qWc.second)
                    {
                        const int iq = std::distance(klist.begin(), std::find(klist.begin(), klist.end(), q_Wc.first));
                        sprintf(fn, "Wcfq_ifreq_%d_iq_%d_I_%zu_J_%zu_id_%d.mtx", ifreq, iq, I, J, mpi_comm_global_h.myid);
                        print_matrix_mm_file(q_Wc.second, Params::output_dir + "/" + fn, 1e-15);
                    }
                }
            }
        }
    }

    LIBRPA::G0W0 s_g0w0(meanfield, kfrac_list, chi0.tfg);
    Profiler::start("g0w0_sigc_IJ", "Build real-space correlation self-energy");
    s_g0w0.build_spacetime(Cs_data, Wc_freq_q, Rlist, period);
    Profiler::stop("g0w0_sigc_IJ");

    // if (Params::debug)
    // { // debug, check sigc_ij
    //     char fn[80];
    //     // ofs_myid << "s_g0w0.sigc_is_f_k_IJ:\n" << s_g0w0.sigc_is_f_k_IJ << "\n";
    //     for (const auto &is_sigc: s_g0w0.sigc_is_f_k_IJ)
    //     {
    //         const int ispin = is_sigc.first;
    //         for (const auto &f_sigc: is_sigc.second)
    //         {
    //             const auto ifreq = chi0.tfg.get_freq_index(f_sigc.first);
    //             for (const auto &k_sigc: f_sigc.second)
    //             {
    //                 const auto &k = k_sigc.first;
    //                 const int ik = std::distance(klist.begin(), std::find(klist.begin(), klist.end(), k));
    //                 for (const auto &I_sigc: k_sigc.second)
    //                 {
    //                     const auto &I = I_sigc.first;
    //                     for (const auto &J_sigc: I_sigc.second)
    //                     {
    //                         const auto &J = J_sigc.first;
    //                         sprintf(fn, "Sigcfq_ispin_%d_ifreq_%d_ik_%d_I_%zu_J_%zu_id_%d.mtx",
    //                                 ispin, ifreq, ik, I, J, mpi_comm_global_h.myid);
    //                         print_matrix_mm_file(J_sigc.second, Params::output_dir + "/" + fn, 1e-15);
    //                     }
    //                 }
    //             }
    //         }
    //     }
    // }

    Profiler::start("g0w0_sigc_rotate_KS", "Construct self-energy in Kohn-Sham space");
    s_g0w0.build_sigc_matrix_KS_kgrid();
    Profiler::stop("g0w0_sigc_rotate_KS");

    if (flag_read_vxc == 0)
    {
        Profiler::start("g0w0_solve_qpe", "Solve quasi-particle equation");
        if (mpi_comm_global_h.is_root())
        {
            std::cout << "Solving quasi-particle equation\n";
        }
        std::vector<cplxdb> imagfreqs;
        for (const auto &freq: chi0.tfg.get_freq_nodes())
        {
            imagfreqs.push_back(cplxdb{0.0, freq});
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
                        for (const auto &freq: chi0.tfg.get_freq_nodes())
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
                                    i_spin+1, i_kpoint+1, i_state+1);
                            e_qp_all[i_spin][i_kpoint][i_state] = std::numeric_limits<double>::quiet_NaN();
                            sigc_all[i_spin][i_kpoint][i_state] = std::numeric_limits<cplxdb>::quiet_NaN();
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
                    printf("spin %2d, k-point %4d: (%.5f, %.5f, %.5f) \n",
                            i_spin+1, i_kpoint+1, k.x, k.y, k.z);
                    printf("%124s\n", banner.c_str());
                    printf("%5s %16s %16s %16s %16s %16s %16s %16s\n", "State", "occ", "e_mf", "v_xc", "v_exx", "ReSigc", "ImSigc", "e_qp");
                    printf("%124s\n", banner.c_str());
                    for (int i_state = 0; i_state < meanfield.get_n_bands(); i_state++)
                    {
                        const auto &occ_state = meanfield.get_weight()[i_spin](i_kpoint, i_state) * meanfield.get_n_kpoints();
                        const auto &eks_state = meanfield.get_eigenvals()[i_spin](i_kpoint, i_state) * HA2EV;
                        const auto &exx_state = exx.Eexx[i_spin][i_kpoint][i_state] * HA2EV;
                        const auto &vxc_state = vxc[i_spin](i_kpoint, i_state) * HA2EV;
                        const auto &resigc = sigc_all[i_spin][i_kpoint][i_state].real() * HA2EV;
                        const auto &imsigc = sigc_all[i_spin][i_kpoint][i_state].imag() * HA2EV;
                        const auto &eqp = e_qp_all[i_spin][i_kpoint][i_state] * HA2EV;
                        printf("%5d %16.5f %16.5f %16.5f %16.5f %16.5f %16.5f %16.5f\n",
                               i_state+1, occ_state, eks_state, vxc_state, exx_state, resigc, imsigc, eqp);
                    }
                    printf("\n");
                }
            }
        }
        Profiler::stop("g0w0_solve_qpe");
    }

    Profiler::start("g0w0_export_sigc_KS", "Export self-energy in KS basis");
    mpi_comm_global_h.barrier();
    if (Params::output_gw_sigc_mat && mpi_comm_global_h.is_root())
    {
        char fn[100];
        for (const auto &ispin_sigc: s_g0w0.sigc_is_ik_f_KS)
        {
            const auto &ispin = ispin_sigc.first;
            for (const auto &ik_sigc: ispin_sigc.second)
            {
                const auto &ik = ik_sigc.first;
                for (const auto &freq_sigc: ik_sigc.second)
                {
                    const auto ifreq = s_g0w0.tfg.get_freq_index(freq_sigc.first);
                    sprintf(fn, "Sigc_fk_mn_ispin_%d_ik_%d_ifreq_%d.mtx", ispin, ik, ifreq);
                    print_matrix_mm_file(freq_sigc.second, Params::output_dir + "/" + fn, 1e-10);
                }
            }
        }
        // for aims analytic continuation reader
        write_self_energy_omega("self_energy_omega.dat", s_g0w0);
    }
    Profiler::stop("g0w0_export_sigc_KS");

    Profiler::stop("g0w0");
}
