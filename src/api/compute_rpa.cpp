// Public API headers
#include "librpa_enums.h"
#include "librpa_func.h"
#include "librpa_compute.h"

// Internal headers
#include "instance_manager.h"
#include "dataset_helper.h"
#include "../core/epsilon.h"
#include "../io/fs.h"
#include "../math/complexmatrix.h"
#include "../utils/profiler.h"
#include "../utils/utils_mem.h"

LIBRPA_C_H_FUNC_WRAP_WOPT(double, librpa_get_rpa_correlation_energy,
                          int n_ibz_kpoints, double *rpa_corr_ibzk_contrib_re,  double *rpa_corr_ibzk_contrib_im)
{
    using namespace librpa_int;
    using librpa_int::global::lib_printf;
    using librpa_int::global::profiler;
    using librpa_int::global::decide_auto_routing;

    double rpa_corr = 0.0;

    profiler.start("api_get_rpa_correlation_energy");
    auto pds = librpa_int::api::get_dataset_instance(h);

    const auto &opts = *p_opts;
    const bool debug = opts.output_level >= LIBRPA_VERBOSE_DEBUG;

    // Prepare time-frequency grids
    initialize_tfgrids(*pds, opts);

    pds->atpairs_local.clear();
    // Determine the actual routing and atom pairs that this process is responsible for
    LibrpaParallelRouting routing = opts.parallel_routing;
    const int n_atoms = pds->atoms.size();
    if (routing == LibrpaParallelRouting::AUTO)
    {
        const auto tot_atpair = generate_atom_pair_from_nat(n_atoms, false);
        routing = decide_auto_routing(tot_atpair.size(), opts.nfreq * pds->pbc.get_n_cells_bvk());
    }
    if(routing == LibrpaParallelRouting::ATOMPAIR || routing == LibrpaParallelRouting::LIBRI)
    {
        auto tri_local_atpair = librpa_int::dispatch_upper_triangular_tasks(
            n_atoms, pds->blacs_ctxt_h.myid, pds->blacs_ctxt_h.nprows, pds->blacs_ctxt_h.npcols,
            pds->blacs_ctxt_h.myprow, pds->blacs_ctxt_h.mypcol);
        for (const auto &p: tri_local_atpair)
            pds->atpairs_local.emplace_back(p);
    }
    else
    {
        pds->atpairs_local = generate_atom_pair_from_nat(n_atoms, false);
    }

    // Initialize response function object
    pds->p_chi0 = std::make_unique<librpa_int::Chi0>(pds->mf, pds->basis_wfc, pds->basis_aux, pds->pbc, pds->tfg, pds->comm_h);
    pds->p_chi0->gf_threshold = opts.gf_threshold;
    pds->p_chi0->libri_threshold_C = opts.libri_chi0_threshold_C;
    pds->p_chi0->libri_threshold_G = opts.libri_chi0_threshold_G;

    auto &chi0 = *(pds->p_chi0);
    // std::cout << "n_abf & " << chi0.atbasis_abf.nb_total;
    // std::cout << "n_abf * " << ds->p_chi0->atbasis_abf.nb_total;

    profiler.start("chi0_build", "Build response function chi0");
    chi0.build(opts.parallel_routing, pds->cs_data, pds->atpairs_local);
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

    const bool use_blacs = opts.use_scalapack_ecrpa && (opts.parallel_routing == LibrpaParallelRouting::ATOMPAIR || opts.parallel_routing == LibrpaParallelRouting::LIBRI);

    if (use_blacs)
    {
        if(pds->mf.get_n_kpoints() == 1)
            corr = compute_RPA_correlation_blacs_2d_gamma_only(chi0, pds->vq, pds->atpairs_local, pds->blacs_ctxt_h);
        else
            corr = compute_RPA_correlation_blacs_2d(chi0, pds->vq, pds->atpairs_local, pds->blacs_ctxt_h);
    }
    else
        corr = compute_RPA_correlation(opts.parallel_routing, *(pds->p_chi0), pds->vq);

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
