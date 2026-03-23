#include "../task.h"
#include "../driver.h"
#include "../read_data.h"
#include "../../src/io/global_io.h"
#include "../../src/io/fs.h"
#include "../../src/io/stl_io_helper.h"
#include "../../src/utils/profiler.h"
#include "../../src/api/instance_manager.h"
#include "../../src/utils/constants.h"
#include "librpa_enums.h"

// #include <fstream>
// #include <sstream>

void driver::task_g0w0_band()
{
    using std::cout;
    using std::endl;
    using std::setw;
    using std::setprecision;
    using std::fixed;
    using std::map;
    using namespace librpa_int;
    using namespace librpa_int::global;

    profiler.start("g0w0_band", "G0W0 quasi-particle band structure calculation");

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
        if (mpi_comm_global_h.is_root())
            std::cout << "Reading dielectric function for head correction" << std::endl;
        std::vector<double> omegas_dielect;
        std::vector<double> dielect_func;
        read_dielec_func(file_df, omegas_dielect, dielect_func);
        ofs_myid << "Dielectric functions read:" << std::endl;
        ofs_myid << "omegas_dielect: " << omegas_dielect << std::endl;
        ofs_myid << "dielect_func:   " << dielect_func << std::endl;
        h.set_dielect_func_imagfreq(omegas_dielect, dielect_func);
    }

    profiler.start("read_vxc", "Load DFT xc potential");
    std::vector<matrix> vxc;
    int flag_read_vxc = read_vxc(driver_params.input_dir + "vxc_out", vxc);
    if (flag_read_vxc == 0)
    {
        if (mpi_comm_global_h.is_root())
            std::cout << "* Success: Read DFT xc potential, solve quasi-particle equation" << std::endl;
    }
    else
        throw LIBRPA_RUNTIME_ERROR("Failed to read DFT xc potential");
    profiler.stop("read_vxc");

    // Build the self-energy matrix (including exchange and correlation)
    h.build_g0w0_sigma(opts);

    const int i_state_low = 0;
    const int i_state_high = n_states;
    const int n_states_calc = i_state_high - i_state_low;
    std::vector<double> vexx_all;
    std::vector<cplxdb> sigc_all;
    {
        const size_t n_local = n_states_calc * n_spins * iks_eigvec_this.size();
        const auto vexx = h.get_exx_pot_kgrid(opts, n_spins, iks_eigvec_this, i_state_low, i_state_high);
        std::vector<double> vxc_flat(n_local);
        ofs_myid << "k_para " << opts.use_kpara_scf_eigvec << endl;
        ofs_myid << "iks_eigvec_this " << iks_eigvec_this << endl;
        ofs_myid << "n_local " << n_local << endl;
        for (int isp = 0; isp != n_spins; isp++)
        {
            const auto start_isp = isp * iks_eigvec_this.size() * n_states_calc;
            for (size_t ik_local = 0; ik_local < iks_eigvec_this.size(); ik_local++)
            {
                const auto ik = iks_eigvec_this[ik_local];
                const auto start_k = start_isp + ik_local * n_states_calc;
                for (int i = 0; i < n_states_calc; i++)
                {
                    vxc_flat[start_k+i] = vxc[isp](ik, i+i_state_low);
                }
            }
        }

        const auto sigc = h.get_g0w0_sigc_kgrid(opts, n_spins, iks_eigvec_this, i_state_low, i_state_high, vxc_flat, vexx);

        if (!opts.use_kpara_scf_eigvec)
        {
            // master process already has all data
            if (myid_global == 0)
            {
                vexx_all = vexx;
                sigc_all = sigc;
            }
        }
        else
        {
            // data parallelized over k-point, collect them to master process
            const size_t n_all = n_states_calc * n_spins * n_kpoints;
            vexx_all.resize(n_all);
            sigc_all.resize(n_all);
            for (int isp = 0; isp != n_spins; isp++)
            {
                const auto st_isp_local = isp * iks_eigvec_this.size() * n_states_calc;
                const auto st_isp = isp * n_kpoints * n_states_calc;
                for (size_t ik_local = 0; ik_local < iks_eigvec_this.size(); ik_local++)
                {
                    const auto ik = iks_eigvec_this[ik_local];
                    const auto st_local = st_isp_local + ik_local * n_states_calc;
                    const auto st = st_isp + ik * n_states_calc;
                    memcpy(sigc_all.data() + st, sigc.data() + st_local, n_states_calc * sizeof(cplxdb));
                    memcpy(vexx_all.data() + st, vexx.data() + st_local, n_states_calc * sizeof(double));
                }
            }
            mpi_comm_global_h.reduce(MPI_IN_PLACE, vexx_all.data(), n_all, 0, MPI_SUM);
            mpi_comm_global_h.reduce(MPI_IN_PLACE, sigc_all.data(), n_all, 0, MPI_SUM);
            if (myid_global != 0)
            {
                vexx_all.clear();
                sigc_all.clear();
            }
        }
    }

    const std::string banner(124, '-');
    const auto &kfrac_list = librpa_int::api::get_dataset_instance(h)->pbc.kfrac_list;
    const auto &mf = librpa_int::api::get_dataset_instance(h)->mf;

    if (flag_read_vxc == 0)
    {
        if (myid_global == 0)
        {
            for (int i_spin = 0; i_spin < n_spins; i_spin++)
            {
                for (int i_kpoint = 0; i_kpoint < n_kpoints; i_kpoint++)
                {
                    const size_t start_k = (i_spin * n_kpoints + i_kpoint) * n_states_calc;
                    for (int i = 0; i < n_states_calc; i++)
                    {
                        const int i_state = i + i_state_low;
                        if (sigc_all[start_k+i].real() == std::numeric_limits<double>::quiet_NaN())
                        {
                            lib_printf("Warning! QPE solver failed for spin %d, kpoint %d, state %d\n",
                                       i_spin+1, i_kpoint+1, i_state+1);
                        }
                    }
                }
            }

            // display results
            lib_printf("Printing quasi-particle energy [unit: eV]\n\n");
            for (int i_spin = 0; i_spin < n_spins; i_spin++)
            {
                for (int i_kpoint = 0; i_kpoint < n_kpoints; i_kpoint++)
                {
                    const auto &k = kfrac_list[i_kpoint];
                    lib_printf("spin %2d, k-point %4d: (%.5f, %.5f, %.5f) \n",
                               i_spin+1, i_kpoint+1, k.x, k.y, k.z);
                    lib_printf("%124s\n", banner.c_str());
                    lib_printf("%5s %16s %16s %16s %16s %16s %16s %16s\n", "State", "occ", "e_mf", "v_xc", "v_exx", "ReSigc", "ImSigc", "e_qp");
                    lib_printf("%124s\n", banner.c_str());
                    const size_t start_k = (i_spin * n_kpoints + i_kpoint) * n_states_calc;
                    for (int i = 0; i < n_states_calc; i++)
                    {
                        const int i_state = i + i_state_low;
                        const auto &occ_state = mf.get_weight()[i_spin](i_kpoint, i_state) * mf.get_n_kpoints();
                        const auto &eks_state = mf.get_eigenvals()[i_spin](i_kpoint, i_state) * HA2EV;
                        const auto &vxc_state = vxc[i_spin](i_kpoint, i_state) * HA2EV;
                        const auto &exx_state = vexx_all[start_k+i] * HA2EV;
                        const auto &resigc = sigc_all[start_k+i].real() * HA2EV;
                        const auto &imsigc = sigc_all[start_k+i].imag() * HA2EV;
                        const auto &eqp = eks_state - vxc_state + exx_state + resigc;
                        lib_printf("%5d %16.5f %16.5f %16.5f %16.5f %16.5f %16.5f %16.5f\n",
                                   i_state+1, occ_state, eks_state, vxc_state, exx_state, resigc, imsigc, eqp);
                    }
                    lib_printf("\n");
                }
            }
        }
    }

    /* Below we handle the band k-points data
     * First load the information of k-points along the k-path */
    profiler.start("g0w0_band_load_band_mf", "Read eigen solutions at band kpoints");
    read_band_kpath_info(driver_params.input_dir + "band_kpath_info");
    const int nkpts_band = kfrac_band.size();
    const int nkpts_band_this = iks_band_eigvec_this.size();
    ofs_myid << "Number of band k-points on this process: " << nkpts_band_this << endl;
    ofs_myid << "iks_band_eigvec_this: " << iks_band_eigvec_this << endl;

    if (mpi_comm_global_h.is_root())
    {
        std::cout << "Band k-points to compute:\n";
        for (int ik = 0; ik < nkpts_band; ik++)
        {
            const auto &k = kfrac_band[ik];
            lib_printf("%5d %12.7f %12.7f %12.7f\n", ik + 1, k.x, k.y, k.z);
        }
    }
    mpi_comm_global_h.barrier();
    read_band_meanfield_data(driver_params.input_dir);
    profiler.stop("g0w0_band_load_band_mf");

    profiler.start("read_vxc_band", "Load DFT xc potential for band");
    const auto vxc_band_all = read_vxc_band(driver_params.input_dir, n_states, n_spins, kfrac_band.size());
    profiler.stop("read_vxc_band");

    // Call APIs to get band EXX potentials and correlation self-energies,
    // then collect them to the master process
    std::vector<double> vexx_band_all;
    std::vector<cplxdb> sigc_band_all;
    {
        const auto &iks_this = iks_band_eigvec_this;
        const size_t n_local = n_states_calc * n_spins * iks_this.size();
        const auto vexx_band = h.get_exx_pot_band_k(opts, n_spins, iks_this, i_state_low, i_state_high);

        std::vector<double> vxc_this(n_local);
        for (int isp = 0; isp != n_spins; isp++)
        {
            const auto start_isp = isp * iks_this.size() * n_states_calc;
            for (size_t ik_this = 0; ik_this < iks_this.size(); ik_this++)
            {
                const auto ik = iks_this[ik_this];
                const auto start_k = start_isp + ik_this * n_states_calc;
                for (int i = 0; i < n_states_calc; i++)
                {
                    vxc_this[start_k+i] = vxc_band_all[isp](ik, i+i_state_low);
                }
            }
        }
        const auto sigc_band = h.get_g0w0_sigc_band_k(opts, n_spins, iks_this, i_state_low, i_state_high, vxc_this, vexx_band);

        profiler.start("collect_exx_sigc_band");
        if (!opts.use_kpara_scf_eigvec)
        {
            // master process already has all data
            if (myid_global == 0)
            {
                vexx_band_all = vexx_band;
                sigc_band_all = sigc_band;
            }
        }
        else
        {
            // data parallelized over k-point, collect them to master process
            const auto n_kpts = n_kpoints_band;
            const size_t n_all = n_states_calc * n_spins * n_kpts;
            vexx_band_all.resize(n_all);
            sigc_band_all.resize(n_all);
            for (int isp = 0; isp != n_spins; isp++)
            {
                const auto st_isp_local = isp * iks_this.size() * n_states_calc;
                const auto st_isp = isp * n_kpts * n_states_calc;
                for (size_t ik_this = 0; ik_this < iks_this.size(); ik_this++)
                {
                    const auto ik = iks_this[ik_this];
                    const auto st_local = st_isp_local + ik_this * n_states_calc;
                    const auto st = st_isp + ik * n_states_calc;
                    memcpy(sigc_band_all.data() + st, sigc_band.data() + st_local, n_states_calc * sizeof(cplxdb));
                    memcpy(vexx_band_all.data() + st, vexx_band.data() + st_local, n_states_calc * sizeof(double));
                }
            }
            mpi_comm_global_h.reduce(MPI_IN_PLACE, vexx_band_all.data(), n_all, 0, MPI_SUM);
            mpi_comm_global_h.reduce(MPI_IN_PLACE, sigc_band_all.data(), n_all, 0, MPI_SUM);
            if (myid_global != 0)
            {
                vexx_band_all.clear();
                sigc_band_all.clear();
            }
        }
        profiler.stop("collect_exx_sigc_band");
    }

    profiler.start("output_g0w0_band");
    const auto &mf_band = librpa_int::api::get_dataset_instance(h)->mf_band;
    {
        const int n_kpts = n_kpoints_band;
        if (myid_global == 0)
        {
            for (int i_spin = 0; i_spin < n_spins; i_spin++)
            {
                for (int i_kpoint = 0; i_kpoint < n_kpts; i_kpoint++)
                {
                    const size_t start_k = (i_spin * n_kpts + i_kpoint) * n_states_calc;
                    for (int i = 0; i < n_states_calc; i++)
                    {
                        const int i_state = i + i_state_low;
                        if (sigc_band_all[start_k+i].real() == std::numeric_limits<double>::quiet_NaN())
                        {
                            lib_printf("Warning! QPE solver failed for spin %d, kpoint %d, state %d\n",
                                       i_spin+1, i_kpoint+1, i_state+1);
                        }
                    }
                }
            }

            for (int i_spin = 0; i_spin < n_spins; i_spin++)
            {
                std::ofstream ofs_ks, ofs_hf, ofs_gw;
                std::stringstream fn_ks, fn_hf, fn_gw;

                fn_gw << "GW_band_spin_" << i_spin + 1 << ".dat";
                fn_hf << "EXX_band_spin_" << i_spin + 1 << ".dat";
                fn_ks << "KS_band_spin_" << i_spin + 1 << ".dat";

                ofs_gw.open(fn_gw.str());
                ofs_hf.open(fn_hf.str());
                ofs_ks.open(fn_ks.str());

                ofs_gw << fixed;
                ofs_hf << fixed;
                ofs_ks << fixed;

                for (int i_kpoint = 0; i_kpoint < n_kpts; i_kpoint++)
                {
                    const auto &k = kfrac_band[i_kpoint];
                    const size_t start_k = (i_spin * n_kpts + i_kpoint) * n_states_calc;

                    ofs_ks << setw(5) << i_kpoint + 1 << setw(15) << setprecision(7) << k.x << setw(15) << setprecision(7) << k.y << std::setw(15) << std::setprecision(7) << k.z;
                    ofs_gw << setw(5) << i_kpoint + 1 << setw(15) << setprecision(7) << k.x << setw(15) << setprecision(7) << k.y << std::setw(15) << std::setprecision(7) << k.z;
                    ofs_hf << setw(5) << i_kpoint + 1 << setw(15) << setprecision(7) << k.x << setw(15) << setprecision(7) << k.y << std::setw(15) << std::setprecision(7) << k.z;
                    for (int i = 0; i < n_states_calc; i++)
                    {
                        const int i_state = i + i_state_low;
                        const auto occ_state = mf_band.get_weight()[i_spin](i_kpoint, i_state) * mf_band.get_n_kpoints();
                        const auto eks_state = mf_band.get_eigenvals()[i_spin](i_kpoint, i_state) * HA2EV;
                        const auto vxc_state = vxc_band_all[i_spin](i_kpoint, i_state) * HA2EV;
                        const auto exx_state = vexx_band_all[start_k+i] * HA2EV;
                        const auto resigc = sigc_band_all[start_k+i].real() * HA2EV;
                        // const auto &imsigc = sigc_band_all[start_k+i].imag() * HA2EV;
                        const auto eqp = eks_state - vxc_state + exx_state + resigc;
                        ofs_ks << setw(15) << std::setprecision(5) << occ_state << setw(15) << setprecision(5) << eks_state;
                        ofs_gw << setw(15) << std::setprecision(5) << occ_state << setw(15) << setprecision(5) << eqp;
                        ofs_hf << setw(15) << std::setprecision(5) << occ_state << setw(15) << setprecision(5) << eks_state - vxc_state + exx_state;
                    }
                    ofs_gw << "\n";
                    ofs_hf << "\n";
                    ofs_ks << "\n";
                }
            }
        }
    }
    profiler.stop("output_g0w0_band");

    // // TODO: parallelize analytic continuation and QPE solver among tasks
    // if (mpi_comm_global_h.is_root())
    // {
    //     const auto &mf = meanfield_band;
    //     const auto n_spin = mf.get_n_spins();
    //     const auto n_kpt = mf.get_n_kpoints();
    //     map<int, map<int, map<int, double>>> e_qp_all;
    //     map<int, map<int, map<int, cplxdb>>> sigc_all;
    //     const auto efermi = mf.get_efermi();
    //
    //     // Spectral function IO handler
    //     std::ofstream wf_sf;
    //     std::vector<cplxdb> omegas;
    //     // FIXME: should check negative n_omegas_sf beforehand for safety and clarity
    //     int n_omegas_sf = (driver_params.sf_omega_end - driver_params.sf_omega_start) / driver_params.sf_omega_step + 1;
    //     const double gf_shift_ha = driver_params.sf_gf_omega_shift / HA2EV;
    //     const double sigc_shift_ha = driver_params.sf_sigc_omega_shift / HA2EV;
    //
    //     if (driver_params.output_gw_spec_func && n_omegas_sf > 0)
    //     {
    //         wf_sf.open("spectral_function_band.dat", std::ios::out | std::ios::binary);
    //
    //         // Initialize the real frequencies
    //         std::vector<double> omegas_db(n_omegas_sf, 0.0);
    //         const double start = driver_params.sf_omega_start;
    //         const double step = driver_params.sf_omega_step;
    //         std::generate(omegas_db.begin(), omegas_db.end(),
    //                       [i = 0, start, step]() mutable
    //                       {
    //                           return start + (i++) * step;
    //                       });
    //         // Output frequency setup
    //         wf_sf.write((char *) &driver_params.sf_omega_start, sizeof(double));
    //         wf_sf.write((char *) &driver_params.sf_omega_end, sizeof(double));
    //         wf_sf.write((char *) &driver_params.sf_omega_step, sizeof(double));
    //         wf_sf.write((char *) &driver_params.sf_gf_omega_shift, sizeof(double));
    //         wf_sf.write((char *) &driver_params.sf_sigc_omega_shift, sizeof(double));
    //         // Output dimension information and Fermi level
    //         wf_sf.write((char *) &n_omegas_sf, sizeof(int));
    //         wf_sf.write((char *) &n_spin, sizeof(int));
    //         wf_sf.write((char *) &n_kpt, sizeof(int));
    //         // Check actual start and end indices of band
    //         int sf_state_start = max(driver_params.sf_state_start, 0);
    //         int sf_state_end = min(driver_params.sf_state_end, mf.get_n_bands() - 1);
    //         wf_sf.write((char *) &sf_state_start, sizeof(int));
    //         wf_sf.write((char *) &sf_state_end, sizeof(int));
    //         wf_sf.write((char *) &efermi, sizeof(double));
    //         // Output the frequencies
    //         wf_sf.write((char *) omegas_db.data(), sizeof(double) * n_omegas_sf);
    //
    //         omegas.resize(n_omegas_sf);
    //         std::transform(omegas_db.cbegin(), omegas_db.cend(), omegas.begin(),
    //                [](double x) { return std::complex<double>(x / HA2EV, 0.0); });
    //     }
    //
    //     for (int i_spin = 0; i_spin < n_spin; i_spin++)
    //     {
    //         for (int i_kpoint = 0; i_kpoint < n_kpt; i_kpoint++)
    //         {
    //             const auto &sigc_sk = s_g0w0.sigc_is_ik_f_KS[i_spin][i_kpoint];
    //             for (int i_state = 0; i_state < mf.get_n_bands(); i_state++)
    //             {
    //                 const auto &eks_state = mf.get_eigenvals()[i_spin](i_kpoint, i_state);
    //                 const auto &exx_state = exx.Eexx[i_spin][i_kpoint][i_state];
    //                 const auto &vxc_state = vxc_band[i_spin](i_kpoint, i_state);
    //                 std::vector<cplxdb> sigc_state;
    //                 for (const auto &freq: chi0.tfg.get_freq_nodes())
    //                 {
    //                     sigc_state.push_back(sigc_sk.at(freq)(i_state, i_state));
    //                 }
    //                 librpa_int::AnalyContPade pade(Params::n_params_anacon, imagfreqs, sigc_state);
    //
    //                 if (driver_params.output_gw_spec_func && n_omegas_sf > 0)
    //                 {
    //                     if (i_state >= driver_params.sf_state_start &&
    //                         i_state <= driver_params.sf_state_end)
    //                     {
    //                         const auto sf = librpa_int::get_specfunc(
    //                             pade, omegas, efermi, eks_state, vxc_state, exx_state,
    //                             sigc_shift_ha, gf_shift_ha);
    //                         wf_sf.write((char *) &eks_state, sizeof(double));
    //                         wf_sf.write((char *) &exx_state, sizeof(double));
    //                         wf_sf.write((char *) &vxc_state, sizeof(double));
    //                         wf_sf.write((char *) sf.data(), sizeof(double) * n_omegas_sf);
    //                     }
    //                 }
    //
    //                 double e_qp;
    //                 cplxdb sigc;
    //                 int flag_qpe_solver = librpa_int::qpe_solver_pade_self_consistent(
    //                     pade, eks_state, efermi, vxc_state, exx_state, e_qp, sigc);
    //                 if (flag_qpe_solver == 0)
    //                 {
    //                     e_qp_all[i_spin][i_kpoint][i_state] = e_qp;
    //                     sigc_all[i_spin][i_kpoint][i_state] = sigc;
    //                 }
    //                 else
    //                 {
    //                     printf("Warning! QPE solver failed for spin %d, kpoint %d, state %d\n",
    //                             i_spin+1, i_kpoint+1, i_state+1);
    //                     e_qp_all[i_spin][i_kpoint][i_state] = std::numeric_limits<double>::quiet_NaN();
    //                     sigc_all[i_spin][i_kpoint][i_state] = std::numeric_limits<cplxdb>::quiet_NaN();
    //                 }
    //             }
    //         }
    //     }
    //
    //     if (driver_params.output_gw_spec_func)
    //     {
    //         wf_sf.close();
    //     }
    //
    //     // display results
    //     for (int i_spin = 0; i_spin < mf.get_n_spins(); i_spin++)
    //     {
    //         std::ofstream ofs_ks;
    //         std::ofstream ofs_hf;
    //         std::ofstream ofs_gw;
    //         std::stringstream fn;
    //
    //         fn << "GW_band_spin_" << i_spin + 1 << ".dat";
    //         ofs_gw.open(fn.str());
    //
    //         fn.str("");
    //         fn.clear();
    //         fn << "EXX_band_spin_" << i_spin + 1 << ".dat";
    //         ofs_hf.open(fn.str());
    //
    //         fn.str("");
    //         fn.clear();
    //         fn << "KS_band_spin_" << i_spin + 1 << ".dat";
    //         ofs_ks.open(fn.str());
    //
    //         ofs_gw << std::fixed;
    //         ofs_hf << std::fixed;
    //         ofs_ks << std::fixed;
    //
    //         for (int i_kpoint = 0; i_kpoint < mf.get_n_kpoints(); i_kpoint++)
    //         {
    //             const auto &k = kfrac_band[i_kpoint];
    //             ofs_ks << std::setw(5) << i_kpoint + 1 << std::setw(15) << std::setprecision(7) << k.x << std::setw(15) << std::setprecision(7) << k.y << std::setw(15) << std::setprecision(7) << k.z;
    //             ofs_gw << std::setw(5) << i_kpoint + 1 << std::setw(15) << std::setprecision(7) << k.x << std::setw(15) << std::setprecision(7) << k.y << std::setw(15) << std::setprecision(7) << k.z;
    //             ofs_hf << std::setw(5) << i_kpoint + 1 << std::setw(15) << std::setprecision(7) << k.x << std::setw(15) << std::setprecision(7) << k.y << std::setw(15) << std::setprecision(7) << k.z;
    //             for (int i_state = 0; i_state < meanfield.get_n_bands(); i_state++)
    //             {
    //                 const auto &occ_state = mf.get_weight()[i_spin](i_kpoint, i_state);
    //                 const auto &eks_state = mf.get_eigenvals()[i_spin](i_kpoint, i_state) * HA2EV;
    //                 const auto &exx_state = exx.Eexx[i_spin][i_kpoint][i_state] * HA2EV;
    //                 const auto &vxc_state = vxc_band[i_spin](i_kpoint, i_state) * HA2EV;
    //                 // const auto &resigc = sigc_all[i_spin][i_kpoint][i_state].real() * HA2EV;
    //                 // const auto &imsigc = sigc_all[i_spin][i_kpoint][i_state].imag() * HA2EV;
    //                 const auto &eqp = e_qp_all[i_spin][i_kpoint][i_state] * HA2EV;
    //                 ofs_ks << std::setw(15) << std::setprecision(5) << occ_state << std::setw(15) << std::setprecision(5) << eks_state;
    //                 ofs_gw << std::setw(15) << std::setprecision(5) << occ_state << std::setw(15) << std::setprecision(5) << eqp;
    //                 ofs_hf << std::setw(15) << std::setprecision(5) << occ_state << std::setw(15) << std::setprecision(5) << eks_state - vxc_state + exx_state;
    //             }
    //             ofs_gw << "\n";
    //             ofs_hf << "\n";
    //             ofs_ks << "\n";
    //         }
    //     }
    // }
    // profiler.stop("g0w0_solve_qpe");


    profiler.stop("g0w0_band");
}
