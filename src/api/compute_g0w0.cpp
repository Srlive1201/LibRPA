// Public API headers
#include "librpa_compute.h"
#include "librpa_enums.h"
#include "librpa_func.h"

// Standard C++ headers
#include <string>

// Headers for API
#include "dataset_helper.h"
#include "instance_manager.h"

// Internal headers
#include "../core/analycont.h"
#include "../core/coulmat.h"
#include "../core/dielecmodel.h"
#include "../core/epsilon.h"
#include "../core/qpe_solver.h"
#include "../io/fs.h"
#include "../math/complexmatrix.h"
#include "../math/utils_matrix_m_mpi.h"
#include "../utils/profiler.h"
#include "../utils/utils_mem.h"

LIBRPA_C_H_FUNC_WRAP_WOPT_NOPAR(void, librpa_build_g0w0_sigma)
{
    using namespace librpa_int;
    using librpa_int::global::profiler;
    using librpa_int::global::lib_printf;

    auto pds = librpa_int::api::get_dataset_instance(h);
    const auto &opts = *p_opts;
    const bool debug = opts.output_level >= LIBRPA_VERBOSE_DEBUG;

    profiler.start("api_build_g0w0_sigma");

    // Prepare time-frequency grids
    initialize_ds_tfgrids(*pds, opts);

    // Decide actual routing
    LibrpaParallelRouting routing = opts.parallel_routing;
    if (routing == LibrpaParallelRouting::AUTO)
    {
        const int n_atoms = pds->atoms.size();
        routing = decide_auto_routing(n_atoms, opts.nfreq * pds->pbc.get_n_cells_bvk());
    }

    // Determine the atom pairs that this process is responsible for
    initialize_ds_atpairs_local(*pds, routing);
    // Redistribute 2D Coulomb matrices to atom-pair blocks if they are parsed
    pds->redistribute_coulomb_blacs2ap();

    // Initialize response function object
    initialize_ds_chi0(*pds, opts);
    auto &chi0 = *(pds->p_chi0);

    profiler.start("chi0_build", "Build response function chi0");
    chi0.build(routing, pds->cs_data, pds->atpairs_local);
    profiler.stop("chi0_build");
    pds->comm_h.barrier();

    if (debug)
    { // debug, check chi0
        for (const auto &chi0q: chi0.get_chi0_q())
        {
            const int ifreq = chi0.tfg.get_freq_index(chi0q.first);
            for (const auto &[q, IJchi0]: chi0q.second)
            {
                const int iq = pds->pbc.get_k_index_full(q);
                for (const auto &[I, J_chi0]: IJchi0)
                {
                    for (const auto &[J, chi0]: J_chi0)
                    {
                        std::stringstream ss;
                        ss << "chi0fq_ifreq_" << ifreq << "_iq_" << iq << "_I_" << I << "_J_" << J << "_id_" << pds->comm_h.myid << ".mtx";
                        print_complex_matrix_mm(chi0, librpa_int::path_as_directory(opts.output_dir) + ss.str(), 1e-15);
                    }
                }
            }
        }
    }
    pds->comm_h.barrier();

    if (!pds->p_exx)
    {
        profiler.start("g0w0_exx", "Build exchange self-energy");
        initialize_ds_exx(*pds, opts);
        profiler.start("ft_vq_cut", "Fourier transform truncated Coulomb");
        const auto VR = librpa_int::FT_Vq(pds->basis_aux, pds->vq_cut, pds->pbc, true);
        profiler.stop("ft_vq_cut");

        profiler.start("g0w0_exx_real_work");
        pds->p_exx->build(routing, pds->basis_aux, pds->cs_data, VR);
        // pds->p_exx->build_KS_kgrid_blacs(pds->blacs_h);
        profiler.stop("g0w0_exx_real_work");
        profiler.stop("g0w0_exx");
        pds->comm_h.barrier();
    }

    profiler.start("g0w0_wc", "Build screened interaction");

    std::vector<double> epsmac_LF_imagfreq_re;
    if (opts.replace_w_head == LIBRPA_SWITCH_ON && pds->epsmacs_imagfreq.size() > 0)
    {
        epsmac_LF_imagfreq_re = interpolate_dielec_func(opts.option_dielect_func, pds->omegas_imagfreq,
                                                        pds->epsmacs_imagfreq, chi0.tfg.get_freq_nodes());
        if (debug)
        {
            if (pds->comm_h.is_root())
            {
                lib_printf("Dielectric function parsed as head correction:\n");
                for (int i = 0; i < opts.nfreq; i++)
                    lib_printf("%d %f %f\n", i+1, chi0.tfg.get_freq_nodes()[i], epsmac_LF_imagfreq_re[i]);
            }
            pds->comm_h.barrier();
        }
    }
    std::vector<std::complex<double>> epsmac_LF_imagfreq(epsmac_LF_imagfreq_re.cbegin(), epsmac_LF_imagfreq_re.cend());

    std::map<double, std::map<Vector3_Order<double>, librpa_int::matrix_m<std::complex<double>>>> Wc_freq_q;
    if (opts.use_scalapack_gw_wc == LIBRPA_SWITCH_ON)
    {
        Wc_freq_q = compute_Wc_freq_q_blacs(chi0, pds->vq, pds->vq_cut, opts.sqrt_coulomb_threshold,
                                            epsmac_LF_imagfreq, pds->blacs_h, pds->desc_abf,
                                            debug, opts.output_dir);
    }
    else
    {
        Wc_freq_q = compute_Wc_freq_q(chi0, pds->vq, pds->vq_cut, opts.sqrt_coulomb_threshold,
                                      epsmac_LF_imagfreq, debug, opts.output_dir);
    }
    profiler.stop("g0w0_wc");

    if (opts.output_level >= LIBRPA_VERBOSE_DEBUG)
    { // debug, check Wc
        for (const auto &[freq, q_Wc]: Wc_freq_q)
        {
            const int ifreq = chi0.tfg.get_freq_index(freq);
            for (const auto &[q, Wc]: q_Wc)
            {
                const int iq = pds->pbc.get_k_index_full(q);
                std::stringstream ss;
                ss << "Wcfq_ifreq_" << ifreq << "_iq_" << iq << ".csc";
                const auto fn = librpa_int::path_as_directory(opts.output_dir) + ss.str();
                write_matrix_elsi_csc_parallel(fn, Wc, pds->desc_abf, 1e-15);
            }
        }
    }

    initialize_ds_g0w0(*pds, opts);
    profiler.start("g0w0_sigc_IJ", "Build real-space correlation self-energy");
    // HACK: choice of space-time is hard-coded. May need to change when more approaches are implemented
    pds->p_g0w0->build_spacetime(routing, pds->basis_aux, pds->cs_data, Wc_freq_q, pds->desc_abf);
    profiler.stop("g0w0_sigc_IJ");

    profiler.stop("api_build_g0w0_sigma");
}

LIBRPA_C_H_FUNC_WRAP_WOPT(void, librpa_get_g0w0_sigc_kgrid,
                          const int n_spins,
                          const int n_kpts_this, const int *iks_this,
                          int i_state_low, int i_state_high, const double *vxc, const double *vexx,
                          double *sigc_re, double *sigc_im)
{
    using namespace librpa_int;
    using librpa_int::global::profiler;
    using librpa_int::global::lib_printf;

    auto pds = librpa_int::api::get_dataset_instance(h);
    const auto &opts = *p_opts;
    const bool debug = opts.output_level >= LIBRPA_VERBOSE_DEBUG;
    i_state_low = std::max(0, i_state_low);
    i_state_high = std::min(pds->mf.get_n_states(), i_state_high);
    if (n_spins != pds->mf.get_n_spins())
    {
        global::ofs_myid << "n_spins != pds->mf.get_n_spins(): " << n_spins << " != " << pds->mf.get_n_spins() << std::endl;
        throw LIBRPA_RUNTIME_ERROR("parsed nspins is not consitent with the SCF starting poing");
    }
    if (i_state_high <= i_state_low) return;
    const int n_states_calc = i_state_high - i_state_low;

    if (!pds->p_g0w0) librpa_build_g0w0_sigma(h, p_opts);

    profiler.start("api_get_g0w0_sigc_kgrid");

    // Decide actual routing
    LibrpaParallelRouting routing = opts.parallel_routing;
    if (routing == LibrpaParallelRouting::AUTO)
    {
        const int n_atoms = pds->atoms.size();
        routing = decide_auto_routing(n_atoms, opts.nfreq * pds->pbc.get_n_cells_bvk());
    }

    profiler.start("g0w0_sigc_rotate_KS", "Correlation self-energy in K-S space");
    pds->p_g0w0->build_sigc_matrix_KS_kgrid_blacs(pds->blacs_h);
    profiler.stop("g0w0_sigc_rotate_KS");

    profiler.start("g0w0_solve_qpe", "Solve quasi-particle equation");
    const auto efermi = pds->mf.get_efermi();
    std::vector<cplxdb> imagfreqs;
    for (const auto &freq: pds->tfg.get_freq_nodes())
    {
        imagfreqs.push_back(cplxdb{0.0, freq});
    }

    for (int isp = 0; isp < n_spins; isp++)
    {
        const int start_isp = isp * n_kpts_this * n_states_calc;
        // ofs_myid << exx_isp << endl;
        for (int ik_this = 0; ik_this < n_kpts_this; ik_this++)
        {
            // When eigenvectors are not parallelized over k, the resulted k-matrices sigc_is_ik_f_KS are collected to master.
	    // Other processes only have dummy data, so we skip AC for them.
            if (!opts.use_kpara_scf_eigvec && pds->blacs_h.myid != 0) continue;

            const int start_k = start_isp + ik_this * n_states_calc;
            const int ik = *(iks_this + ik_this);
            global::ofs_myid << "Start QPE solver for spin " << isp + 1
                             << " kpoint " << ik + 1 << std::endl;
            for (int i = 0; i < n_states_calc; i++)
            {
                const int i_state = i + i_state_low;
                const auto &eks_state = pds->mf.get_eigenvals()[isp](ik, i_state);
                const auto &exx_state = vexx[start_k+i];
                const auto &vxc_state = vxc[start_k+i];
                std::vector<cplxdb> sigc_state;
                const auto &sigc_isp = pds->p_g0w0->sigc_is_ik_f_KS[isp];
                auto it = sigc_isp.find(ik);
                if (it == sigc_isp.cend())
                    throw LIBRPA_RUNTIME_ERROR("fail to locate sigc at ik = " + std::to_string(ik));
                for (const auto &[freq, mat]: it->second)
                {
                    sigc_state.emplace_back(mat(i_state, i_state));
                }
                double e_qp;
                cplxdb sigc;
                sigc_re[start_k+i] = std::numeric_limits<double>::quiet_NaN();
                sigc_im[start_k+i] = std::numeric_limits<double>::quiet_NaN();
                librpa_int::AnalyContPade pade(opts.n_params_anacon, imagfreqs, sigc_state);
                int flag_qpe_solver = librpa_int::qpe_solver_pade_self_consistent(
                    pade, eks_state, efermi, vxc_state, exx_state, e_qp, sigc);
                if (flag_qpe_solver == 0)
                {
                    sigc_re[start_k+i] = sigc.real();
                    sigc_im[start_k+i] = sigc.imag();
                }
                else
                {
                    global::ofs_myid << "Warning! QPE solver failed for spin " << isp + 1
                                     << " kpoint " << ik + 1 << " state " << i_state + 1
                                     << std::endl;
                }
            }
        }
    }
    profiler.stop("g0w0_solve_qpe");

    profiler.stop("api_get_g0w0_sigc_kgrid");
}

LIBRPA_C_H_FUNC_WRAP_WOPT(void, librpa_get_g0w0_sigc_band_k,
                          const int n_spins,
                          const int n_kpts_band_this, const int *iks_band_this,
                          int i_state_low, int i_state_high, const double *vxc_band, const double *vexx_band,
                          double *sigc_band_re, double *sigc_band_im)
{
    using namespace librpa_int;
    using librpa_int::global::profiler;
    using librpa_int::global::lib_printf;

    auto pds = librpa_int::api::get_dataset_instance(h);
    const auto &opts = *p_opts;
    const bool debug = opts.output_level >= LIBRPA_VERBOSE_DEBUG;
    i_state_low = std::max(0, i_state_low);
    i_state_high = std::min(pds->mf_band.get_n_states(), i_state_high);
    const int n_spins_band = pds->mf_band.get_n_spins();
    if (n_spins != n_spins_band)
    {
        global::ofs_myid << "n_spins != pds->mf_band.get_n_spins(): " << n_spins << " != " << n_spins_band << std::endl;
        throw LIBRPA_RUNTIME_ERROR("parsed n_spins is not consitent with the band input data");
    }
    if (i_state_high <= i_state_low) return;
    const int n_states_calc = i_state_high - i_state_low;

    if (!pds->p_g0w0) librpa_build_g0w0_sigma(h, p_opts);

    profiler.start("api_get_g0w0_sigc_band_k");

    // Decide actual routing
    LibrpaParallelRouting routing = opts.parallel_routing;
    if (routing == LibrpaParallelRouting::AUTO)
    {
        const int n_atoms = pds->atoms.size();
        routing = decide_auto_routing(n_atoms, opts.nfreq * pds->pbc.get_n_cells_bvk());
    }

    profiler.start("g0w0_sigc_rotate_KS", "Correlation self-energy in K-S space");
    pds->p_g0w0->build_sigc_matrix_KS_band_blacs(pds->mf_band.get_eigenvectors(),
                                                 pds->kfrac_band_list, pds->atoms, pds->blacs_h);
    profiler.stop("g0w0_sigc_rotate_KS");

    profiler.start("g0w0_solve_qpe", "Solve quasi-particle equation");
    const auto efermi = pds->mf_band.get_efermi();
    std::vector<cplxdb> imagfreqs;
    for (const auto &freq: pds->tfg.get_freq_nodes())
    {
        imagfreqs.push_back(cplxdb{0.0, freq});
    }

    for (int isp = 0; isp < n_spins; isp++)
    {
        const int start_isp = isp * n_kpts_band_this * n_states_calc;
        // ofs_myid << exx_isp << endl;
        for (int ik_this = 0; ik_this < n_kpts_band_this; ik_this++)
        {
            // When eigenvectors are not parallelized over k, the resulted k-matrices sigc_is_ik_f_KS are collected to master.
	    // Other processes only have dummy data, so we skip AC for them.
            if (!opts.use_kpara_scf_eigvec && pds->blacs_h.myid != 0) continue;

            const int start_k = start_isp + ik_this * n_states_calc;
            const int ik = *(iks_band_this + ik_this);
            global::ofs_myid << "Start QPE solver for spin " << isp + 1
                             << " kpoint " << ik + 1 << std::endl;
            for (int i = 0; i < n_states_calc; i++)
            {
                const int i_state = i + i_state_low;
                const auto &eks_state = pds->mf_band.get_eigenvals()[isp](ik, i_state);
                const auto &exx_state = vexx_band[start_k+i];
                const auto &vxc_state = vxc_band[start_k+i];
                std::vector<cplxdb> sigc_state;
                const auto &sigc_isp = pds->p_g0w0->sigc_is_ik_f_KS[isp];
                auto it = sigc_isp.find(ik);
                if (it == sigc_isp.cend())
                    throw LIBRPA_RUNTIME_ERROR("fail to locate sigc at ik = " + std::to_string(ik));
                for (const auto &[freq, mat]: it->second)
                {
                    sigc_state.emplace_back(mat(i_state, i_state));
                }
                double e_qp;
                cplxdb sigc;
                sigc_band_re[start_k+i] = std::numeric_limits<double>::quiet_NaN();
                sigc_band_im[start_k+i] = std::numeric_limits<double>::quiet_NaN();
                librpa_int::AnalyContPade pade(opts.n_params_anacon, imagfreqs, sigc_state);
                int flag_qpe_solver = librpa_int::qpe_solver_pade_self_consistent(
                    pade, eks_state, efermi, vxc_state, exx_state, e_qp, sigc);
                if (flag_qpe_solver == 0)
                {
                    sigc_band_re[start_k+i] = sigc.real();
                    sigc_band_im[start_k+i] = sigc.imag();
                }
                else
                {
                    global::ofs_myid << "Warning! QPE solver failed for spin " << isp + 1
                                     << " kpoint " << ik + 1 << " state " << i_state + 1
                                     << std::endl;
                }
            }
        }
    }
    profiler.stop("g0w0_solve_qpe");

    profiler.stop("api_get_g0w0_sigc_band_k");
}
