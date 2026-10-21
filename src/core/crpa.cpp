#include "crpa.h"

#include "../io/global_io.h"
#include "../utils/error.h"
#include "../utils/profiler.h"
#include "crpa_workflow.h"

namespace librpa_int
{

CrpaResult compute_crpa_onsite(CrpaContext &context, const Chi0::BandSelection &selection,
                               const std::vector<CorrelatedSubspace *> &frames)
{
    global::profiler.start("crpa_onsite", "Onsite constrained and full-RPA interactions");
    try
    {
        const int unpaired_optical =
            (context.full_headwing == nullptr) != (context.constrained_headwing == nullptr);
        int any_unpaired_optical = 0;
        context.chi0.comm_h.allreduce(&unpaired_optical, &any_unpaired_optical, 1, MPI_MAX);
        if (any_unpaired_optical)
            throw LIBRPA_INVALID_ARGUMENT(
                "standalone cRPA requires matching full and constrained optical models");
        auto projection = prepare_crpa_local_projection(context, selection, frames);
        auto full_response = build_crpa_total_response(context);
        auto excluded = build_crpa_excluded_response(context, selection);
        auto constrained =
            form_crpa_constrained_response(context, full_response, std::move(excluded));
        global::lib_printf_root("cRPA: screening Pr=P0-Pd for U\n");
        {
            auto wc =
                screen_crpa_response(context, std::move(constrained), context.constrained_headwing);
            auto values = contract_crpa_interaction(context, projection, wc);
            for (std::size_t tensor = 0; tensor < projection.bare_result.tensors.size(); ++tensor)
                projection.bare_result.tensors[tensor].u = std::move(values[tensor]);
        }
        // Release constrained Wc before constructing the full-RPA kernels.
        global::lib_printf_root("cRPA: screening P0 for full-RPA W\n");
        auto wc = screen_crpa_response(context, std::move(full_response), context.full_headwing);
        auto values = contract_crpa_interaction(context, projection, wc);
        auto result = std::move(projection.bare_result);
        for (std::size_t tensor = 0; tensor < result.tensors.size(); ++tensor)
            result.tensors[tensor].w = std::move(values[tensor]);
        global::profiler.stop("crpa_onsite");
        return result;
    }
    catch (...)
    {
        global::profiler.stop("crpa_onsite");
        throw;
    }
}

}  // namespace librpa_int
