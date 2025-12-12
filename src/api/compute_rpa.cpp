// Public API headers
#include "librpa_enums.h"
#include "librpa_func.h"
#include "librpa_compute.h"

// Internal headers
#include "../core/epsilon.h"
#include "../io/fs.h"
// #include "../io/stl_io_helper.h"
#include "../math/complexmatrix.h"
#include "../utils/profiler.h"
#include "../utils/utils_mem.h"
#include "dataset_helper.h"
#include "instance_manager.h"

LIBRPA_C_H_FUNC_WRAP_WOPT(void, librpa_get_imaginary_frequency_grids,
                          double *omegas, double *weights)
{
    using namespace librpa_int;
    using librpa_int::global::lib_printf;
    using librpa_int::global::profiler;

    profiler.start("api_get_imaginary_freq_grids");

    auto pds = librpa_int::api::get_dataset_instance(h);
    const auto &opts = *p_opts;
    initialize_ds_tfgrids(*pds, opts);
    const auto v_freqs = pds->tfg.get_freq_nodes();
    const auto v_wegihts = pds->tfg.get_freq_weights();
    // global::ofs_myid << v_freqs << std::endl;
    memcpy(omegas, v_freqs.data(), opts.nfreq * sizeof(double));
    memcpy(weights, v_wegihts.data(), opts.nfreq * sizeof(double));

    profiler.stop("api_get_imaginary_freq_grids");
}

LIBRPA_C_H_FUNC_WRAP_WOPT(double, librpa_get_rpa_correlation_energy,
                          int n_ibz_kpoints, double *rpa_corr_ibzk_contrib_re,  double *rpa_corr_ibzk_contrib_im)
{
    using namespace librpa_int;
    using librpa_int::global::lib_printf;
    using librpa_int::global::profiler;

    double rpa_corr = 0.0;

    profiler.start("api_get_rpa_correlation_energy");
    auto pds = librpa_int::api::get_dataset_instance(h);
    const auto &opts = *p_opts;

    const auto ks_local = pds->mf.get_iks_local();
    const bool seems_k_para = as_int(ks_local.size()) < pds->mf.get_n_kpoints();
    if (seems_k_para)
    {
        global::ofs_myid << "KS eigenvectors of SCF start point seem distributed over k-points" << std::endl;
        if (opts.use_kpara_scf_eigvec == LIBRPA_SWITCH_ON)
        {
            global::ofs_myid << "Option use_kpara_scf_eigvec is on: input consistent" << std::endl;
        }
        else
        {
            global::ofs_myid << "Option use_kpara_scf_eigvec is off while the eigenvectors are distributed" << std::endl;
            throw LIBRPA_RUNTIME_ERROR("inconsistent input");
        }
    }

    const bool debug = opts.output_level >= LIBRPA_VERBOSE_DEBUG;

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

    // std::cout << "n_abf & " << chi0.atbasis_abf.nb_total;
    // std::cout << "n_abf * " << pds->p_chi0->atbasis_abf.nb_total;

    profiler.start("chi0_build", "Build response function chi0");
    chi0.build(routing, pds->cs_data, pds->atpairs_local);
    profiler.stop("chi0_build");

    if (debug)
    { // debug, check chi0
        profiler.start("chi0_write_matrix_mm", "Export chi0 matrices");
        char fn[80];
        for (const auto &[freq, q_IJchi]: chi0.get_chi0_q())
        {
            const int ifreq = chi0.tfg.get_freq_index(freq);
            for (const auto &[q, I_Jchi]: q_IJchi)
            {
                const int iq = pds->pbc.get_k_index_ibz(q);
                for (const auto &[I, J_chi]: I_Jchi)
                {
                    for (const auto &[J, chi]: J_chi)
                    {
                        sprintf(fn, "chi0fq_ifreq_%d_iq_%d_I_%zu_J_%zu_id_%d.mtx", ifreq, iq, I, J, pds->comm_h.myid);
                        print_complex_matrix_mm(chi, path_as_directory(opts.output_dir) + fn);
                    }
                }
            }
        }
        profiler.stop("chi0_write_matrix_mm");
    }

    profiler.start("chi0_release_free");
    release_free_mem();
    pds->comm_h.barrier();
    profiler.stop("chi0_release_free");

    profiler.start("EcRPA", "Compute RPA correlation Energy");
    CorrEnergy corr;

    const bool use_blacs = opts.use_scalapack_ecrpa && (routing == LibrpaParallelRouting::ATOMPAIR || routing == LibrpaParallelRouting::LIBRI);

    if (use_blacs)
    {
        if(pds->mf.get_n_kpoints() == 1)
            corr = compute_RPA_correlation_blacs_2d_gamma_only(chi0, pds->vq, pds->atpairs_local, pds->blacs_h);
        else
            corr = compute_RPA_correlation_blacs_2d(chi0, pds->vq, pds->atpairs_local, pds->blacs_h);
    }
    else
        corr = compute_RPA_correlation(routing, *(pds->p_chi0), pds->vq);

    rpa_corr = corr.value.real();
    if (corr.qcontrib.size() != as_size(n_ibz_kpoints))
    {
        if (pds->comm_h.is_root())
        {
            lib_printf(
                "WARNING: parsed n_ibz_kpoints is not consistent with generated CorrEne object: %d != %zu\n",
                as_size(n_ibz_kpoints), corr.qcontrib.size());
        }
    }

    for (const auto &[q, corr_q]: corr.qcontrib)
    {
        const auto iq = pds->pbc.get_k_index_ibz(q);
        if (iq == -1)
            throw LIBRPA_RUNTIME_ERROR("Internal error, failed to find irreducible k-point");
        rpa_corr_ibzk_contrib_re[iq] = corr_q.real();
        rpa_corr_ibzk_contrib_im[iq] = corr_q.imag();
    }
    profiler.stop("EcRPA");

    // Works done, free the object and reset the pointer to nullptr
    pds->p_chi0.reset();
    profiler.stop("api_get_rpa_correlation_energy");

    return rpa_corr;
}
