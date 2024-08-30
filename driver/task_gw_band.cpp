#include "task_gw_band.h"

#include <sstream>
#include <fstream>

#include "envs_mpi.h"
#include "envs_io.h"
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
#include "read_data.h"
#include "write_aims.h"
#include "driver_utils.h"

void task_g0w0_band()
{
    using LIBRPA::envs::mpi_comm_global_h;
    using LIBRPA::envs::ofs_myid;
    using LIBRPA::utils::lib_printf;

    Profiler::start("g0w0_band", "G0W0 quasi-particle band structure calculation");

    Vector3_Order<int> period {kv_nmp[0], kv_nmp[1], kv_nmp[2]};
    auto Rlist = construct_R_grid(period);

    vector<Vector3_Order<double>> qlist;
    for ( auto q_weight: irk_weight)
    {
        qlist.push_back(q_weight.first);
    }

    Chi0 chi0(meanfield, klist, Params::nfreq);
    chi0.gf_R_threshold = Params::gf_R_threshold;

    Profiler::start("chi0_build", "Build response function chi0");
    chi0.build(Cs_data, Rlist, period, local_atpair, qlist,
               TFGrids::get_grid_type(Params::tfgrids_type), true);
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

    std::vector<double> epsmac_LF_imagfreq_re;

    if (Params::replace_w_head)
    {
        std::vector<double> omegas_dielect;
        std::vector<double> dielect_func;
        read_dielec_func("dielecfunc_out", omegas_dielect, dielect_func);

        epsmac_LF_imagfreq_re = interpolate_dielec_func(
                Params::option_dielect_func, omegas_dielect, dielect_func,
                chi0.tfg.get_freq_nodes());
    }

    Profiler::start("g0w0_exx", "Build exchange self-energy");
    auto exx = LIBRPA::Exx(meanfield, kfrac_list);
    {
        Profiler::start("ft_vq_cut", "Fourier transform truncated Coulomb");
        const auto VR = FT_Vq(Vq_cut, Rlist, true);
        Profiler::stop("ft_vq_cut");

        Profiler::start("g0w0_exx_real_work");
        exx.build(Cs_data, Rlist, period, VR);
        Profiler::stop("g0w0_exx_real_work");
    }
    Profiler::stop("g0w0_exx");
    std::flush(ofs_myid);

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
    std::flush(ofs_myid);

    /* Below we handle the band k-points data
     * First load the information of k-points along the k-path */
    Profiler::start("g0w0_band_load_band_mf", "Read eigen solutions at band kpoints");
    int n_basis_band, n_states_band, n_spin_band;
    std::vector<Vector3_Order<double>> kfrac_band =
        read_band_kpath_info(n_basis_band, n_states_band, n_spin_band);
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

    auto meanfield_band = read_meanfield_band(
            n_basis_band, n_states_band, n_spin_band, kfrac_band.size());

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
    auto vxc_band = read_vxc_band(n_states_band, n_spin_band, kfrac_band.size());
    Profiler::stop("read_vxc");
    std::flush(ofs_myid);

    std::vector<cplxdb> imagfreqs;
    for (const auto &freq: chi0.tfg.get_freq_nodes())
    {
        imagfreqs.push_back(cplxdb{0.0, freq});
    }

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
        for (int i_spin = 0; i_spin < mf.get_n_spins(); i_spin++)
        {
            std::ofstream ofs_hf;
            std::ofstream ofs_gw;
            std::stringstream fn;
            fn << "GW_band_spin_" << i_spin + 1 << ".dat";
            ofs_gw.open(fn.str());
            fn.str("");
            fn.clear();

            fn << "EXX_band_spin_" << i_spin + 1 << ".dat";
            ofs_hf.open(fn.str());

            ofs_gw << std::fixed;
            ofs_hf << std::fixed;

            for (int i_kpoint = 0; i_kpoint < mf.get_n_kpoints(); i_kpoint++)
            {
                const auto &k = kfrac_band[i_kpoint];
                ofs_gw << std::setw(5) << i_kpoint + 1 << std::setw(15) << std::setprecision(7) << k.x << std::setw(15) << std::setprecision(7) << k.y << std::setw(15) << std::setprecision(7) << k.z;
                ofs_hf << std::setw(5) << i_kpoint + 1 << std::setw(15) << std::setprecision(7) << k.x << std::setw(15) << std::setprecision(7) << k.y << std::setw(15) << std::setprecision(7) << k.z;
                for (int i_state = 0; i_state < meanfield.get_n_bands(); i_state++)
                {
                    const auto &occ_state = mf.get_weight()[i_spin](i_kpoint, i_state);
                    const auto &eks_state = mf.get_eigenvals()[i_spin](i_kpoint, i_state) * HA2EV;
                    const auto &exx_state = exx.Eexx[i_spin][i_kpoint][i_state] * HA2EV;
                    const auto &vxc_state = vxc_band[i_spin](i_kpoint, i_state) * HA2EV;
                    // const auto &resigc = sigc_all[i_spin][i_kpoint][i_state].real() * HA2EV;
                    // const auto &imsigc = sigc_all[i_spin][i_kpoint][i_state].imag() * HA2EV;
                    const auto &eqp = e_qp_all[i_spin][i_kpoint][i_state] * HA2EV;
                    ofs_gw << std::setw(15) << std::setprecision(5) << occ_state << std::setw(15) << std::setprecision(5) << eqp;
                    ofs_hf << std::setw(15) << std::setprecision(5) << occ_state << std::setw(15) << std::setprecision(5) << eks_state - vxc_state + exx_state;
                }
                ofs_gw << "\n";
                ofs_hf << "\n";
            }
        }
    }
    Profiler::stop("g0w0_solve_qpe");


    Profiler::stop("g0w0_band");
}
