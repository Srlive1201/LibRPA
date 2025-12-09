// #include "../src/core/exx.h"
// #include "../src/core/gw.h"
// #include "../src/core/meanfield.h"
// #include "../src/core/ri.h"
// #include "../src/core/utils_timefreq.h"
// #include "../src/math/utils_matrix_m_mpi.h"
// #include "../src/mpi/envs_blacs.h"
// #include "../src/mpi/global_mpi.h"
// #include "../src/utils/constants.h"

#include <memory>
#include <sstream>
#include "../task.h"
#include "../driver.h"
#include "../read_data.h"
#include "../write_aims.h"
#include "../../src/core/analycont.h"
#include "../../src/core/qpe_solver.h"
#include "../../src/core/epsilon.h"
#include "../../src/io/global_io.h"
#include "../../src/io/fs.h"
#include "../../src/core/dielecmodel.h"
#include "../../src/utils/profiler.h"
#include "../../src/core/coulmat.h"
#include "../../src/api/instance_manager.h"
#include "../../src/utils/constants.h"
#include "../../src/math/utils_matrix_m_mpi.h"

#include "librpa_enums.h"
// #include "driver_params.h"
// #include "driver_utils.h"
// #include "read_data.h"
// #include "write_aims.h"

void driver::task_g0w0()
{
    using std::cout;
    using std::endl;
    using std::map;
    using namespace librpa_int::global;
    using librpa_int::global::lib_printf;
    using librpa_int::cplxdb;
    using librpa_int::HA2EV;

    profiler.start("g0w0", "G0W0 quasi-particle calculation");

    profiler.start("read_vq_cut", "Load truncated Coulomb");

    // Prepare input specific to G0W0
    auto routing = opts.parallel_routing;
    if (routing == LibrpaParallelRouting::AUTO) routing = librpa_int::decide_auto_routing(n_atoms, n_kpoints * opts.nfreq);

    if (routing == LibrpaParallelRouting::RTAU)
    {
        read_Vq_full(driver_params.input_dir, "coulomb_cut_", true);
    }
    else
    {
        // NOTE: local_atpair set during read_data::read_ri.
        read_Vq_row(driver_params.input_dir, "coulomb_cut_", opts.vq_threshold, local_atpair, true);
    }
    profiler.stop("read_vq_cut");

    const auto file_df = driver_params.input_dir + "dielecfunc_out";
    if (driver::get_bool(opts.replace_w_head) && librpa_int::path_exists(file_df.c_str()))
    {
        std::vector<double> omegas_dielect;
        std::vector<double> dielect_func;
        read_dielec_func(file_df, omegas_dielect, dielect_func);
        h.set_dielect_func_imagfreq(omegas_dielect, dielect_func);
    }

    // Build the self-energy matrix
    h.build_g0w0_sigma(opts);

    auto pds = librpa_int::api::get_dataset_instance(h.get_c_handler());
    auto &chi0 = *(pds->p_chi0);

    // if (opts.output_level >= LIBRPA_VERBOSE_DEBUG)
    // { // debug, check chi0
    //     for (const auto &chi0q: chi0.get_chi0_q())
    //     {
    //         const int ifreq = chi0.tfg.get_freq_index(chi0q.first);
    //         for (const auto &[q, IJchi0]: chi0q.second)
    //         {
    //             const int iq = ds->pbc.get_k_index_full(q);
    //             for (const auto &[I, J_chi0]: IJchi0)
    //             {
    //                 for (const auto &[J, chi0]: J_chi0)
    //                 {
    //                     std::stringstream ss;
    //                     ss << "chi0fq_ifreq_" << ifreq << "_iq_" << iq << "_I_" << I << "_J_" << J << "_id_" << ds->comm_h.myid << ".mtx";
    //                     print_complex_matrix_mm(chi0, librpa_int::path_as_directory(opts.output_dir) + ss.str(), 1e-15);
    //                 }
    //             }
    //         }
    //     }
    // }

    profiler.start("read_vxc", "Load DFT xc potential");
    std::vector<matrix> vxc;
    int flag_read_vxc = read_vxc(driver_params.input_dir + "vxc_out", vxc);
    if (mpi_comm_global_h.myid == 0)
    {
        if (flag_read_vxc == 0)
        {
            std::cout << "* Success: Read DFT xc potential, will solve quasi-particle equation\n";
        }
        else
        {
            std::cout << "*   Error: Read DFT xc potential, switch off solving quasi-particle equation\n";
        }
    }
    profiler.stop("read_vxc");
    std::flush(ofs_myid);

    // // Since EXX contribution will be printed along with sigma_c when Vxc is read
    // // we print it here when Vxc read failed
    // if (flag_read_vxc != 0 && mpi_comm_global_h.is_root())
    // {
    //     for (int isp = 0; isp != ds->mf.get_n_spins(); isp++)
    //     {
    //         lib_printf("Spin channel %1d\n", isp+1);
    //         for (int ik = 0; ik != ds->mf.get_n_kpoints(); ik++)
    //         {
    //             cout << "k-point " << ik + 1 << ": " << ds->pbc.klist[ik] << endl;
    //             lib_printf("%-4s  %-10s  %-10s\n", "Band", "e_exx (Ha)", "e_exx (eV)");
    //             for (int ib = 0; ib != ds->mf.get_n_bands(); ib++)
    //                 lib_printf("%4d  %10.5f  %10.5f\n", ib+1, ds->p_exx->Eexx[isp][ik][ib],
    //                            HA2EV * ds->p_exx->Eexx[isp][ik][ib]);
    //             lib_printf("\n");
    //         }
    //     }
    // }
    // mpi_comm_global_h.barrier();

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

    profiler.start("g0w0_sigc_rotate_KS", "Construct self-energy in Kohn-Sham space");
    pds->p_g0w0->build_sigc_matrix_KS_kgrid_blacs(pds->blacs_ctxt_h);
    profiler.stop("g0w0_sigc_rotate_KS");

    if (flag_read_vxc == 0)
    {
        profiler.start("g0w0_solve_qpe", "Solve quasi-particle equation");
        if (pds->comm_h.is_root())
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
            const auto efermi = pds->mf.get_efermi();
            for (int i_spin = 0; i_spin < pds->mf.get_n_spins(); i_spin++)
            {
                for (int i_kpoint = 0; i_kpoint < pds->mf.get_n_kpoints(); i_kpoint++)
                {
                    const auto &sigc_sk = pds->p_g0w0->sigc_is_ik_f_KS[i_spin][i_kpoint];
                    // ofs_myid << sigc_sk.size() << endl;
                    for (int i_state = 0; i_state < pds->mf.get_n_bands(); i_state++)
                    {
                        const auto &eks_state = pds->mf.get_eigenvals()[i_spin](i_kpoint, i_state);
                        const auto &exx_state = pds->p_exx->Eexx[i_spin][i_kpoint][i_state];
                        const auto &vxc_state = vxc[i_spin](i_kpoint, i_state);
                        std::vector<cplxdb> sigc_state;
                        for (const auto &freq: chi0.tfg.get_freq_nodes())
                        {
                            sigc_state.push_back(sigc_sk.at(freq)(i_state, i_state));
                        }
                        librpa_int::AnalyContPade pade(opts.n_params_anacon, imagfreqs, sigc_state);
                        double e_qp;
                        cplxdb sigc;
                        int flag_qpe_solver = librpa_int::qpe_solver_pade_self_consistent(
                            pade, eks_state, efermi, vxc_state, exx_state, e_qp, sigc);
                        if (flag_qpe_solver == 0)
                        {
                            e_qp_all[i_spin][i_kpoint][i_state] = e_qp;
                            sigc_all[i_spin][i_kpoint][i_state] = sigc;
                        }
                        else
                        {
                            lib_printf("Warning! QPE solver failed for spin %d, kpoint %d, state %d\n",
                                       i_spin+1, i_kpoint+1, i_state+1);
                            e_qp_all[i_spin][i_kpoint][i_state] = std::numeric_limits<double>::quiet_NaN();
                            sigc_all[i_spin][i_kpoint][i_state] = std::numeric_limits<cplxdb>::quiet_NaN();
                        }
                    }
                }
            }

            // display results
            const std::string banner(124, '-');
            lib_printf("Printing quasi-particle energy [unit: eV]\n\n");
            for (int i_spin = 0; i_spin < pds->mf.get_n_spins(); i_spin++)
            {
                for (int i_kpoint = 0; i_kpoint < pds->mf.get_n_kpoints(); i_kpoint++)
                {
                    const auto &k = pds->pbc.kfrac_list[i_kpoint];
                    lib_printf("spin %2d, k-point %4d: (%.5f, %.5f, %.5f) \n",
                               i_spin+1, i_kpoint+1, k.x, k.y, k.z);
                    lib_printf("%124s\n", banner.c_str());
                    lib_printf("%5s %16s %16s %16s %16s %16s %16s %16s\n", "State", "occ", "e_mf", "v_xc", "v_exx", "ReSigc", "ImSigc", "e_qp");
                    lib_printf("%124s\n", banner.c_str());
                    for (int i_state = 0; i_state < pds->mf.get_n_bands(); i_state++)
                    {
                        const auto &occ_state = pds->mf.get_weight()[i_spin](i_kpoint, i_state) * pds->mf.get_n_kpoints();
                        const auto &eks_state = pds->mf.get_eigenvals()[i_spin](i_kpoint, i_state) * HA2EV;
                        const auto &exx_state = pds->p_exx->Eexx[i_spin][i_kpoint][i_state] * HA2EV;
                        const auto &vxc_state = vxc[i_spin](i_kpoint, i_state) * HA2EV;
                        const auto &resigc = sigc_all[i_spin][i_kpoint][i_state].real() * HA2EV;
                        const auto &imsigc = sigc_all[i_spin][i_kpoint][i_state].imag() * HA2EV;
                        const auto &eqp = e_qp_all[i_spin][i_kpoint][i_state] * HA2EV;
                        lib_printf("%5d %16.5f %16.5f %16.5f %16.5f %16.5f %16.5f %16.5f\n",
                                   i_state+1, occ_state, eks_state, vxc_state, exx_state, resigc, imsigc, eqp);
                    }
                    lib_printf("\n");
                }
            }
        }
        profiler.stop("g0w0_solve_qpe");
    }

    profiler.start("g0w0_export_sigc_KS", "Export self-energy in KS basis");
    pds->comm_h.barrier();
    if (driver::get_bool(opts.output_gw_sigc_mat) && pds->comm_h.is_root())
    {
        char fn[100];
        for (const auto &ispin_sigc: pds->p_g0w0->sigc_is_ik_f_KS)
        {
            const auto &ispin = ispin_sigc.first;
            for (const auto &ik_sigc: ispin_sigc.second)
            {
                const auto &ik = ik_sigc.first;
                for (const auto &freq_sigc: ik_sigc.second)
                {
                    const auto ifreq = pds->p_g0w0->tfg.get_freq_index(freq_sigc.first);
                    sprintf(fn, "Sigc_fk_mn_ispin_%d_ik_%d_ifreq_%d.mtx", ispin, ik, ifreq);
                    print_matrix_mm_file(freq_sigc.second, librpa_int::path_as_directory(opts.output_dir) + fn, 1e-10);
                }
            }
        }
        // for aims analytic continuation reader
    }
    // Limit writing of file to the master process (it has all the data)
    if (pds->comm_h.is_root())
    {
        write_self_energy_omega("self_energy_omega.dat", *(pds->p_g0w0));
    }
    profiler.stop("g0w0_export_sigc_KS");

    profiler.stop("g0w0");
}
