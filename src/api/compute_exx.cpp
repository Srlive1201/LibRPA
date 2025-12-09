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

LIBRPA_C_H_FUNC_WRAP_WOPT_NOPAR(void, librpa_build_exx)
{
    using namespace librpa_int;
    using librpa_int::global::profiler;
    using librpa_int::global::lib_printf;

    auto pds = librpa_int::api::get_dataset_instance(h);
    const auto &opts = *p_opts;
    const bool debug = opts.output_level >= LIBRPA_VERBOSE_DEBUG;

    profiler.start("api_build_exx");

    // Decide actual routing
    LibrpaParallelRouting routing = opts.parallel_routing;
    if (routing == LibrpaParallelRouting::AUTO)
    {
        const int n_atoms = pds->atoms.size();
        routing = decide_auto_routing(n_atoms, opts.nfreq * pds->pbc.get_n_cells_bvk());
    }

    initialize_ds_exx(*pds, opts);
    {
        profiler.start("ft_vq_cut", "Fourier transform truncated Coulomb");
        const auto VR = librpa_int::FT_Vq(pds->basis_aux, pds->vq_cut, pds->pbc, true);
        profiler.stop("ft_vq_cut");

        profiler.start("exx_real_work");
        pds->p_exx->build(routing, pds->basis_aux, pds->cs_data, VR);
        // pds->p_exx->build_KS_kgrid_blacs(pds->blacs_ctxt_h);
        profiler.stop("exx_real_work");
    }

    profiler.stop("api_build_exx");
}
