// Public API headers
#include "librpa_compute.h"

// Internal headers
#include "../core/coulmat.h"
#include "../io/fs.h"
#include "../io/global_io.h"
#include "../io/stl_io_helper.h"
#include "../utils/error.h"
#include "../utils/profiler.h"
#include "../utils/utils_mem.h"
#include "compute_helper.h"
#include "dataset_helper.h"
#include "instance_manager.h"

#include <algorithm>
#include <map>
#include <string>
#include <vector>

namespace
{

// Map a compact collected k-list back to its buffer position.
std::map<int, int> make_ik_pos_map(const std::vector<int> &iks)
{
    std::map<int, int> ik_pos;
    for (int i = 0; i != static_cast<int>(iks.size()); ++i)
    {
        ik_pos.emplace(iks[i], i);
    }
    return ik_pos;
}

// Publish EXX diagonal values from the rank that owns each rotated k-block and
// fill the caller's output buffer in its original requested k-point order.
void collect_exx_diag_to_callers(const std::map<int, std::map<int, std::map<int, double>>> &Eexx,
                                 const librpa_int::MpiCommHandler &comm_h,
                                 const bool publish_local_values,
                                 const int n_spins,
                                 const int n_kpoints,
                                 const int n_kpts_this,
                                 const int *iks_this,
                                 const int i_state_low,
                                 const int n_states_calc,
                                 double *vexx)
{
    const auto iks_collect =
        librpa_int::api::collect_requested_iks(comm_h, n_kpts_this, iks_this, n_kpoints);
    if (iks_collect.empty()) return;

    const auto ik_pos = make_ik_pos_map(iks_collect);
    const int nk_collect = static_cast<int>(iks_collect.size());
    const int n_owner = n_spins * nk_collect;
    const int n_value = n_owner * n_states_calc;

    std::vector<int> owners(n_owner, 0);
    std::vector<int> owners_sum(n_owner, 0);
    std::vector<double> values(n_value, 0.0);
    std::vector<double> values_sum(n_value, 0.0);

    if (publish_local_values)
    {
        for (int isp = 0; isp != n_spins; ++isp)
        {
            const auto it_sp = Eexx.find(isp);
            if (it_sp == Eexx.cend()) continue;
            for (int ik_collect = 0; ik_collect != nk_collect; ++ik_collect)
            {
                const int ik = iks_collect[ik_collect];
                const auto it_k = it_sp->second.find(ik);
                if (it_k == it_sp->second.cend()) continue;

                owners[isp * nk_collect + ik_collect] = 1;
                for (int i = 0; i != n_states_calc; ++i)
                {
                    const int idx = (isp * nk_collect + ik_collect) * n_states_calc + i;
                    values[idx] = it_k->second.at(i_state_low + i);
                }
            }
        }
    }

    comm_h.allreduce(values.data(), values_sum.data(), n_value, MPI_SUM);
    comm_h.allreduce(owners.data(), owners_sum.data(), n_owner, MPI_SUM);

    for (int isp = 0; isp != n_spins; ++isp)
    {
        const int start_isp = isp * n_kpts_this * n_states_calc;
        for (int ik_this = 0; ik_this != n_kpts_this; ++ik_this)
        {
            const int ik = iks_this[ik_this];
            const int ik_collect = ik_pos.at(ik);
            const int owner_count = owners_sum[isp * nk_collect + ik_collect];
            if (owner_count != 1)
            {
                throw LIBRPA_RUNTIME_ERROR(
                    "failed to locate a unique EXX value owner for spin = " +
                    std::to_string(isp) + " ik = " + std::to_string(ik) +
                    ", owner count = " + std::to_string(owner_count));
            }
            const int start_k = start_isp + ik_this * n_states_calc;
            for (int i = 0; i != n_states_calc; ++i)
            {
                const int idx = (isp * nk_collect + ik_collect) * n_states_calc + i;
                vexx[start_k + i] = values_sum[idx];
            }
        }
    }
}

} // namespace

void librpa_build_exx(LibrpaHandler* h, const LibrpaOptions *p_opts)
{
    using namespace librpa_int;
    using librpa_int::global::profiler;
    using librpa_int::global::lib_printf;

    auto pds = librpa_int::api::get_dataset_instance(h);
    const auto &opts = *p_opts;
    // const bool debug = opts.output_level >= LIBRPA_VERBOSE_DEBUG;
    pds->is_band_calc_done = false;
    initialize_ds_global_ddla(*pds, opts);

    profiler.start("api_build_exx");

    // Decide actual routing
    LibrpaParallelRouting routing = opts.parallel_routing;
    if (routing == LIBRPA_ROUTING_AUTO)
    {
        const int n_atoms = pds->atoms.size();
        routing = decide_auto_routing(n_atoms, opts.nfreq * pds->pbc.get_n_cells_bvk());
    }

    if (opts.use_kpara_scf_eigvec == LIBRPA_SWITCH_ON)
    {
        if (routing == LIBRPA_ROUTING_LIBRI)
            pds->redistribute_eigvecs_kpara_2d();
        else
            pds->redistribute_eigvecs_kpara();
    }

    // Determine the atom pairs that this process is responsible for
    initialize_ds_atpairs_local(*pds, routing);
    // Redistribute 2D Coulomb matrices to atom-pair blocks if they are parsed
    pds->redistribute_coulomb_blacs2ap();

    initialize_ds_exx(*pds, opts);
    const bool use_shrink_abfs = opts.use_shrink_abfs == LIBRPA_SWITCH_ON;
    const auto &basis_aux_exx = use_shrink_abfs ? pds->basis_aux_shrink : pds->basis_aux;
    const auto &cs_data_exx = use_shrink_abfs ? pds->cs_data_shrink : pds->cs_data;
    const auto &coul = opts.use_fullcoul_exx ? pds->vq : pds->vq_cut;
    profiler.start("ft_vq_cut", "Fourier transform truncated Coulomb");
    bool coul_mat_is_replicated = false;
    const auto VR = librpa_int::FT_Vq(
        pds->comm_h, basis_aux_exx, pds->symmetry_context, coul, pds->pbc, true,
        opts.use_symmetry_exx == LIBRPA_SWITCH_ON, &coul_mat_is_replicated);
    profiler.stop("ft_vq_cut");

    profiler.start("exx_real_work");
    pds->p_exx->build(routing, basis_aux_exx, cs_data_exx, VR, coul_mat_is_replicated);
    // pds->p_exx->build_KS_kgrid_blacs(pds->blacs_h);
    profiler.stop("exx_real_work");
    // global::ofs_myid << pds->p_exx->exx_IJR << std::endl;
    pds->comm_h.barrier();

    profiler.stop("api_build_exx");
}

void librpa_get_exx_pot_kgrid(LibrpaHandler *h, const LibrpaOptions *p_opts, const int n_spins,
                              const int n_kpts_this, const int *iks_this, int i_state_low,
                              int i_state_high, double *vexx)
{
    using std::endl;
    using namespace librpa_int;
    using librpa_int::global::profiler;
    using librpa_int::global::ofs_myid;
    using librpa_int::global::lib_printf;

    auto pds = librpa_int::api::get_dataset_instance(h);
    const auto &opts = *p_opts; // TODO: add a flag to control whether to use blacs or lapack
    initialize_ds_global_ddla(*pds, opts);
    i_state_low = std::max(0, i_state_low);
    i_state_high = std::min(pds->mf.get_n_states(), i_state_high);
    if (n_spins != pds->mf.get_n_spins())
    {
        global::ofs_myid << "n_spins != pds->mf.get_n_spins(): " << n_spins << " != " << pds->mf.get_n_spins() << endl;
        throw LIBRPA_RUNTIME_ERROR("parsed nspins is not consitent with the SCF starting point");
    }
    if (i_state_high <= i_state_low)
    {
        return;
    }

    if (!pds->p_exx)
    {
        librpa_build_exx(h, p_opts);
    }

    // const bool debug = opts.output_level >= LIBRPA_VERBOSE_DEBUG;

    profiler.start("api_get_exx_pot_kgrid");
    auto &pexx = pds->p_exx;
    // ofs_myid << pexx->exx_IJR << endl;
    // TODO: make choosing blacs/non-blacs method a run time option
    pexx->build_KS_kgrid_blacs(pds->blacs_h, opts.use_gpu_replace_scalapack);
    pds->is_band_calc_done = false;
    if (opts.output_exx_ks_mat_k == LIBRPA_SWITCH_ON)
    {
        profiler.start("exx_export_KS", "Export exact exchange in KS basis");
        pexx->write_exx_matrices_KS_binary(
            opts.output_dir, "kgrid",
            opts.istate_output_mat_start, opts.istate_output_mat_end);
        profiler.stop("exx_export_KS");
    }
    const int n_states_calc = i_state_high - i_state_low;
    const bool publish_local_values =
        opts.use_kpara_scf_eigvec == LIBRPA_SWITCH_ON || pds->blacs_h.myid == 0;
    collect_exx_diag_to_callers(pexx->Eexx, pds->comm_h, publish_local_values,
                                n_spins, pds->mf.get_n_kpoints(), n_kpts_this, iks_this,
                                i_state_low, n_states_calc, vexx);
    profiler.stop("api_get_exx_pot_kgrid");
}

void librpa_get_exx_pot_band_k(LibrpaHandler *h, const LibrpaOptions *p_opts, const int n_spins,
                               const int n_kpts_band_this, const int *iks_band_this,
                               int i_state_low, int i_state_high, double *vexx_band)
{
    using std::endl;
    using namespace librpa_int;
    using librpa_int::global::profiler;
    using librpa_int::global::ofs_myid;
    using librpa_int::global::lib_printf;

    auto pds = librpa_int::api::get_dataset_instance(h);
    const auto &opts = *p_opts; // TODO: add a flag to control whether to use blacs or lapack
    initialize_ds_global_ddla(*pds, opts);

    if (!pds->is_band_data_set || pds->mf_band.get_n_spins() == 0)
        throw LIBRPA_RUNTIME_ERROR("Meanfield data for band calculation is not set");

    i_state_low = std::max(0, i_state_low);
    i_state_high = std::min(pds->mf.get_n_states(), i_state_high);
    if (n_spins != pds->mf.get_n_spins())
    {
        global::ofs_myid << "n_spins != pds->mf_band.get_n_spins(): " << n_spins << " != " << pds->mf_band.get_n_spins() << endl;
        throw LIBRPA_RUNTIME_ERROR("parsed nspins is not consitent with the SCF starting point");
    }
    if (i_state_high <= i_state_low)
    {
        return;
    }

    if (!pds->p_exx)
    {
        librpa_build_exx(h, p_opts);
    }

    // const bool debug = opts.output_level >= LIBRPA_VERBOSE_DEBUG;

    profiler.start("api_get_exx_pot_band_k");
    auto &pexx = pds->p_exx;
    // ofs_myid << pexx->exx_IJR << endl;
    // TODO: make choosing blacs/non-blacs method a run time option
    if (!pds->is_band_calc_done || pexx->Eexx.empty())
    {
        if (!pds->is_band_calc_done && pds->p_g0w0)
        {
            pds->p_g0w0->reset_kspace();
        }
        pexx->reset_kspace();
        const auto bvk_remap = librpa_int::api::build_band_bvk_remap(
            pds->atoms, pds->pbc, opts.option_bvk_remap);
        if (opts.use_kpara_scf_eigvec == LIBRPA_SWITCH_ON)
        {
            pds->redistribute_band_eigvecs_kpara_2d();
        }
        pexx->build_KS_band_blacs(pds->mf_band.get_eigenvectors(), pds->kfrac_band_list,
                                  bvk_remap, pds->blacs_h, opts.use_gpu_replace_scalapack,
                                  pds->band_data_id);
        pds->is_band_calc_done = true;
    }
    if (opts.output_exx_ks_mat_k == LIBRPA_SWITCH_ON)
    {
        profiler.start("exx_export_KS", "Export exact exchange in KS basis");
        pexx->write_exx_matrices_KS_binary(
            opts.output_dir, "band_" + std::to_string(pds->band_data_id),
            opts.istate_output_mat_start, opts.istate_output_mat_end);
        profiler.stop("exx_export_KS");
    }
    const int n_states_calc = i_state_high - i_state_low;
    const bool publish_local_values =
        opts.use_kpara_scf_eigvec == LIBRPA_SWITCH_ON || pds->blacs_h.myid == 0;
    collect_exx_diag_to_callers(pexx->Eexx, pds->comm_h, publish_local_values,
                                n_spins, pds->mf_band.get_n_kpoints(), n_kpts_band_this,
                                iks_band_this, i_state_low, n_states_calc, vexx_band);
    profiler.stop("api_get_exx_pot_band_k");
}
