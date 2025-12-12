// Public API headers
#include "librpa_enums.h"
#include "librpa_func.h"
#include "librpa_compute.h"

// Internal headers
#include "../core/coulmat.h"
#include "../core/dielecmodel.h"
#include "../core/epsilon.h"
#include "../io/fs.h"
#include "../math/complexmatrix.h"
#include "../math/utils_matrix_m_mpi.h"
#include "../utils/profiler.h"
#include "../utils/utils_mem.h"
#include "dataset_helper.h"
#include "instance_manager.h"

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

    profiler.start("g0w0_exx", "Build exchange self-energy");
    initialize_ds_exx(*pds, opts);
    {
        profiler.start("ft_vq_cut", "Fourier transform truncated Coulomb");
        const auto VR = librpa_int::FT_Vq(pds->basis_aux, pds->vq_cut, pds->pbc, true);
        profiler.stop("ft_vq_cut");

        profiler.start("g0w0_exx_real_work");
        pds->p_exx->build(routing, pds->basis_aux, pds->cs_data, VR);
        // pds->p_exx->build_KS_kgrid_blacs(pds->blacs_h);
        profiler.stop("g0w0_exx_real_work");
    }
    profiler.stop("g0w0_exx");

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
