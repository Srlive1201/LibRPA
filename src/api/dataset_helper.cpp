// Public API headers
#include "librpa_enums.h"

// Internal headers
#include "../utils/profiler.h"
#include "dataset_helper.h"

namespace librpa_int
{

void initialize_ds_tfgrids(Dataset &ds, const LibrpaOptions &opts)
{
    global::profiler.start("initialize_ds_tfgrids");
    ds.tfg.reset(opts.nfreq);
    double emin = opts.tfgrids_freq_min;
    double eintv = opts.tfgrids_freq_interval;
    double emax = opts.tfgrids_freq_max;
    double tmin = opts.tfgrids_time_min;
    double tintv = opts.tfgrids_time_interval;
    if (opts.tfgrids_type == LibrpaTimeFreqGrid::Minimax)
    {
        ds.mf.get_E_min_max(emin, emax);
    }
    ds.tfg.generate(opts.tfgrids_type, emin, eintv, emax, tmin, tintv);
    global::profiler.stop("initialize_ds_tfgrids");
}

void initialize_ds_atpairs_local(Dataset &ds, LibrpaParallelRouting routing)
{
    const int n_atoms = ds.atoms.size();

    ds.atpairs_local.clear();

    if (routing == LibrpaParallelRouting::AUTO)
    {
        throw LIBRPA_RUNTIME_ERROR("internal error: routing should be decided before initialize_ds_atpairs_local, not AUTO");
    }
    else if(routing == LibrpaParallelRouting::ATOMPAIR || routing == LibrpaParallelRouting::LIBRI)
    {
        auto tri_local_atpair = librpa_int::dispatch_upper_triangular_tasks(
            n_atoms, ds.blacs_h.myid, ds.blacs_h.nprows, ds.blacs_h.npcols,
            ds.blacs_h.myprow, ds.blacs_h.mypcol);
        for (const auto &p: tri_local_atpair)
            ds.atpairs_local.emplace_back(p);
    }
    else
    {
        ds.atpairs_local = generate_atom_pair_from_nat(n_atoms, false);
    }
}

void initialize_ds_exx(Dataset &ds, const LibrpaOptions &opts) noexcept
{
    global::profiler.start("initialize_ds_exx");
    const bool is_eigvec_k_distributed = opts.use_kpara_scf_eigvec == LIBRPA_SWITCH_ON;
    ds.p_exx = std::make_unique<librpa_int::Exx>(ds.mf, ds.basis_wfc, ds.pbc, ds.comm_h,
                                                 is_eigvec_k_distributed);
    ds.p_exx->libri_threshold_C = opts.libri_exx_threshold_C;
    ds.p_exx->libri_threshold_D = opts.libri_exx_threshold_D;
    ds.p_exx->libri_threshold_V = opts.libri_exx_threshold_V;
    global::profiler.stop("initialize_ds_exx");
}

void initialize_ds_chi0(Dataset &ds, const LibrpaOptions &opts) noexcept
{
    global::profiler.start("initialize_ds_chi0");
    const bool is_eigvec_k_distributed = opts.use_kpara_scf_eigvec == LIBRPA_SWITCH_ON;
    ds.p_chi0 = std::make_unique<librpa_int::Chi0>(ds.mf, ds.basis_wfc, ds.basis_aux, ds.pbc,
                                                   ds.tfg, ds.comm_h, is_eigvec_k_distributed);
    ds.p_chi0->gf_threshold = opts.gf_threshold;
    ds.p_chi0->libri_threshold_C = opts.libri_chi0_threshold_C;
    ds.p_chi0->libri_threshold_G = opts.libri_chi0_threshold_G;
    global::profiler.stop("initialize_ds_chi0");
}

void initialize_ds_g0w0(Dataset &ds, const LibrpaOptions &opts) noexcept
{
    global::profiler.start("initialize_ds_g0w0");
    const bool is_eigvec_k_distributed = opts.use_kpara_scf_eigvec == LIBRPA_SWITCH_ON;
    ds.p_g0w0 = std::make_unique<librpa_int::G0W0>(ds.mf, ds.basis_wfc, ds.pbc, ds.tfg, ds.comm_h,
                                                   is_eigvec_k_distributed);
    ds.p_g0w0->libri_threshold_C = opts.libri_g0w0_threshold_C;
    ds.p_g0w0->libri_threshold_G = opts.libri_g0w0_threshold_G;
    ds.p_g0w0->libri_threshold_Wc = opts.libri_g0w0_threshold_Wc;
    ds.p_g0w0->output_dir = opts.output_dir;
    ds.p_g0w0->output_sigc_mat = opts.output_gw_sigc_mat == LIBRPA_SWITCH_ON;
    ds.p_g0w0->output_sigc_mat_rt = opts.output_gw_sigc_mat_rt == LIBRPA_SWITCH_ON;
    ds.p_g0w0->output_sigc_mat_rf = opts.output_gw_sigc_mat_rf == LIBRPA_SWITCH_ON;
    global::profiler.stop("initialize_ds_g0w0");
}

}
