#pragma once

#include <memory>

#include "crpa.h"
#include "crpa_window_input.h"

namespace librpa_int
{

//! Owning model for one electronic state; rebuild after changing its E/C/f.
struct CrpaPreparedModel
{
    CrpaWindowBuildResult windows;
    std::vector<std::unique_ptr<CorrelatedSubspace>> frames;
    Chi0::BandSelection selection;
    bool retain_raw_model = true;
};

//! Local orbital-pair vertices in the active auxiliary space and retained bare v.
struct CrpaLocalProjection
{
    CrpaResult bare_result;
    std::vector<std::vector<std::map<Vector3_Order<double>, ComplexMatrix>>> vertices;
};

using CrpaScreenedResponse = std::map<double, std::map<Vector3_Order<double>, Matz>>;

CrpaPreparedModel prepare_crpa_model(const AtomicBasis &basis, const BasisConvention &convention,
                                     const std::vector<int> &atom_types,
                                     const std::vector<std::string> &species_labels,
                                     const std::vector<ComplexMatrix> &spatial_overlap,
                                     const PeriodicBoundaryData &pbc, const MeanField &meanfield,
                                     const CrpaWindowOptions &options);

CrpaLocalProjection prepare_crpa_local_projection(CrpaContext &context,
                                                  const Chi0::BandSelection &selection,
                                                  const std::vector<CorrelatedSubspace *> &frames);

//! Response stages preserve the resident Chi0 response and its incoming selector.
Chi0QMap build_crpa_total_response(CrpaContext &context);
Chi0QMap build_crpa_excluded_response(CrpaContext &context, const Chi0::BandSelection &selection);
Chi0QMap form_crpa_constrained_response(CrpaContext &context, const Chi0QMap &full_response,
                                        Chi0QMap excluded_response);

//! Consume an owning response and return Wc in the active 2D-distributed layout.
//! Optical option 3 and magnetic full-q restoration retain the standalone conventions.
CrpaScreenedResponse screen_crpa_response(CrpaContext &context, Chi0QMap response,
                                          diele_func *optical = nullptr);

//! Return retained bare v + Wc, indexed [tensor][frequency], without changing Wc.
std::vector<std::vector<ComplexMatrix>> contract_crpa_interaction(
    CrpaContext &context, const CrpaLocalProjection &projection,
    const CrpaScreenedResponse &screened_response);

/*! Compute U using the caller's P0, without building P0 or screening full W.
 * The caller must supply P0 for the same current E/C/f, minimax nodes, q view,
 * active auxiliary basis and MPI atom-pair ownership represented by context.
 * The borrowed P0 (including context.chi0.get_chi0_q()), resident response map
 * and incoming band selector are preserved. Transient Green-function workspaces
 * may be rebuilt. Returned tensors contain bare and U; their W vectors are empty.
 */
CrpaResult evaluate_crpa_onsite_from_response(CrpaContext &context,
                                              const Chi0::BandSelection &selection,
                                              const std::vector<CorrelatedSubspace *> &frames,
                                              const Chi0QMap &full_response);

}  // namespace librpa_int
