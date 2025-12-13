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

static void collect_atpairs_all(Dataset &ds)
{
    const auto &comm_h = ds.comm_h;
    const auto &atpairs_local = ds.atpairs_local;
    const int np_this = atpairs_local.size();
    std::vector<int> np_all(comm_h.nprocs, 0);
    np_all[comm_h.myid] = np_this;
    int np_max;
    comm_h.allreduce(&np_this, &np_max, 1, MPI_MAX);
    comm_h.allreduce(MPI_IN_PLACE, np_all.data(), comm_h.nprocs, MPI_SUM);
    if (np_max == 0) return;
    std::vector<size_t> pairs_all(comm_h.nprocs * np_max * 2, 0);
    const int st = comm_h.myid * np_max * 2;
    for (int ip = 0; ip < np_this; ip++)
    {
        pairs_all[st + ip * 2] = atpairs_local[ip].first;
        pairs_all[st + ip * 2 + 1] = atpairs_local[ip].second;
    }
    comm_h.allreduce(MPI_IN_PLACE, pairs_all.data(), pairs_all.size(), MPI_SUM);
    for (int pid = 0; pid < comm_h.nprocs; pid++)
    {
        if (pid == comm_h.myid)
        {
            std::set<atpair_t> atpairs_this(atpairs_local.cbegin(), atpairs_local.cend());
            ds.atpairs_unique_all.emplace(pid, atpairs_this);
        }
        else
        {
            const int np_this = np_all[pid];
            std::set<atpair_t> atpairs_this;
            const int st = pid * np_max * 2;
            for (int ip = 0; ip < np_this; ip++)
            {
                atpairs_this.insert({pairs_all[st + ip * 2], pairs_all[st + ip * 2 + 1]});
            }
            ds.atpairs_unique_all.emplace(pid, atpairs_this);
        }
    }
}

void initialize_ds_atpairs_local(Dataset &ds, LibrpaParallelRouting routing)
{
    global::profiler.start(__FUNCTION__);

    ds.atpairs_local.clear();
    const int n_atoms_basis_wfc = ds.basis_wfc.n_atoms;
    const int n_atoms_basis_aux = ds.basis_aux.n_atoms;
    const int n_atoms_struc = ds.atoms.size();
    const int n_atoms = n_atoms_struc > 0? n_atoms_struc : std::max(n_atoms_basis_aux, n_atoms_basis_wfc);
    if (n_atoms == 0)
        throw LIBRPA_RUNTIME_ERROR("Number of atoms can not be extracted, please set structure or basis first");

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
    collect_atpairs_all(ds);
    global::profiler.stop(__FUNCTION__);
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
