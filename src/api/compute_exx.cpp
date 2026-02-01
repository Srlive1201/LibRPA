// Public API headers
#include "librpa_enums.h"
#include "librpa_func.h"
#include "librpa_compute.h"

// Internal headers
#include "../core/coulmat.h"
#include "../io/fs.h"
#include "../io/stl_io_helper.h"
#include "../utils/error.h"
#include "../utils/profiler.h"
#include "../utils/utils_mem.h"
#include "dataset_helper.h"
#include "instance_manager.h"

LIBRPA_C_H_FUNC_WRAP_WOPT_NOPAR(void, librpa_build_exx)
{
    using namespace librpa_int;
    using librpa_int::global::profiler;
    using librpa_int::global::lib_printf;

    auto pds = librpa_int::api::get_dataset_instance(h);
    const auto &opts = *p_opts;
    // const bool debug = opts.output_level >= LIBRPA_VERBOSE_DEBUG;

    profiler.start("api_build_exx");

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

    initialize_ds_exx(*pds, opts);
    profiler.start("ft_vq_cut", "Fourier transform truncated Coulomb");
    const auto VR = librpa_int::FT_Vq(pds->basis_aux, pds->vq_cut, pds->pbc, true);
    profiler.stop("ft_vq_cut");

    profiler.start("exx_real_work");
    pds->p_exx->build(routing, pds->basis_aux, pds->cs_data, VR);
    // pds->p_exx->build_KS_kgrid_blacs(pds->blacs_h);
    profiler.stop("exx_real_work");
    // global::ofs_myid << pds->p_exx->exx_IJR << std::endl;
    pds->comm_h.barrier();

    profiler.stop("api_build_exx");
}

LIBRPA_C_H_FUNC_WRAP_WOPT(void, librpa_get_exx_pot_kgrid, const int n_spins, const int n_kpoints_local,
                          const int *iks_local, int i_state_low, int i_state_high, double *vexx)
{
    using std::endl;
    using namespace librpa_int;
    using librpa_int::global::profiler;
    using librpa_int::global::ofs_myid;
    using librpa_int::global::lib_printf;

    auto pds = librpa_int::api::get_dataset_instance(h);
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

    // const auto &opts = *p_opts; // TODO: add a flag to control whether to use blacs or lapack
    // const bool debug = opts.output_level >= LIBRPA_VERBOSE_DEBUG;

    profiler.start("api_get_exx_pot_kgrid");
    auto &pexx = pds->p_exx;
    // ofs_myid << pexx->exx_IJR << endl;
    // TODO: make choosing blacs/non-blacs method a run time option
    pexx->build_KS_kgrid_blacs(pds->blacs_h);
    const int n_states_calc = i_state_high - i_state_low;
    for (int isp = 0; isp < n_spins; isp++)
    {
        const auto it = pexx->Eexx.find(isp);
        // NOTE: Without use_kpara_scf_eigvec, exx data is only collected to master.
        // This policy might need to be adjusted later.
        if (it == pexx->Eexx.cend()) continue;
        const int start_isp = isp * n_kpoints_local * n_states_calc;
        // ofs_myid << exx_isp << endl;
        for (int ik_local = 0; ik_local < n_kpoints_local; ik_local++)
        {
            const int start_k = start_isp + ik_local * n_states_calc;
            const int ik = iks_local[ik_local];
            const auto it_k = it->second.find(ik);
            if (it_k == it->second.cend()) continue;
            const auto &exx_k = it_k->second;
            for (int i = 0; i < n_states_calc; i++)
            {
                vexx[start_k+i] = exx_k.at(i_state_low+i);
            }
        }
    }
    profiler.stop("api_get_exx_pot_kgrid");
}
