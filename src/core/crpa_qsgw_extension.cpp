#include "crpa_qsgw_extension.h"

#include <algorithm>
#include <optional>
#include <utility>

#include "../io/global_io.h"
#include "../mpi/collective_stage.h"
#include "../utils/error.h"

namespace librpa_int
{
CrpaQsgwExtension::CrpaQsgwExtension(CrpaTaskInput input, CrpaQsgwControls controls,
                                     CrpaQsgwOutputSink sink)
    : input_(std::move(input)), controls_(std::move(controls)), sink_(std::move(sink))
{
    if (controls_.interval <= 0 || !sink_)
        throw LIBRPA_INVALID_ARGUMENT("QSGW cRPA requires a positive interval and output sink");
    if (controls_.track_orbitals && controls_.tracking_definition.empty())
        throw LIBRPA_INVALID_ARGUMENT("QSGW cRPA tracking requires a model definition");
}

void CrpaQsgwExtension::initialize(Dataset&, const LibrpaOptions& options, int)
{
    if (options.read_sigc_mat_rf == LIBRPA_SWITCH_ON)
        throw LIBRPA_INVALID_ARGUMENT(
            "QSGW cRPA requires the live response; read_sigc_mat_rf skips that build");
    options_ = options;
}

bool CrpaQsgwExtension::wants_response(int entering_state) const
{
    return entering_state % controls_.interval == 0;
}

void CrpaQsgwExtension::update_tracking(const CrpaTaskResult& result)
{
    if (!controls_.track_orbitals) return;
    tracking_ = snapshot_crpa_tracking(result.model);
    input_.window.previous_parent_spin_k = tracking_.parent_phi;
}

void CrpaQsgwExtension::emit(Dataset& dataset, CrpaTaskResult result, int state, int gw_update,
                             bool final)
{
    update_tracking(result);
    std::optional<OrbitalWeightResult> analysis;
    if (controls_.output_analysis)
        collective_all_stage(dataset.comm_h, "cRPA orbital analysis",
                             [&]
                             {
                                 std::vector<int> atom_types(dataset.basis_wfc.n_atoms);
                                 for (std::size_t atom = 0; atom < atom_types.size(); ++atom)
                                     atom_types[atom] = dataset.atoms.types.at(atom);
                                 analysis = analyze_orbital_weights(
                                     dataset.basis_wfc, dataset.basis_convention, atom_types,
                                     dataset.mf, input_.spatial_overlap, dataset.pbc.weight_k,
                                     input_.window.residual_tol);
                             });
    sink_(result, analysis ? &*analysis : nullptr, state, gw_update, final);
}

void CrpaQsgwExtension::evaluate_response(Dataset& dataset, const LibrpaOptions& options,
                                          const Chi0QMap& response, int entering_state,
                                          int gw_update)
{
    auto model = prepare_crpa_task_model(dataset, options, input_, CrpaTaskMode::SharedGw);
    emit(dataset, evaluate_crpa_task_from_response(dataset, options, model, response),
         entering_state, gw_update, false);
}

void CrpaQsgwExtension::state_updated(Dataset& dataset, int)
{
    if (!controls_.track_orbitals) return;
    update_tracking(prepare_crpa_task(dataset, options_, input_));
    const double capture =
        *std::min_element(tracking_.parent_capture.begin(), tracking_.parent_capture.end());
    global::lib_printf_root(
        "cRPA parent transport: minimum retained probability %.16e; no response build\n", capture);
}

void CrpaQsgwExtension::read_checkpoint(Dataset& dataset, const std::string& filename,
                                        int iteration)
{
    if (!controls_.track_orbitals) return;
    collective_all_stage(dataset.comm_h, "cRPA parent tracking restart",
                         [&]
                         {
                             tracking_ = read_crpa_tracking_checkpoint(
                                 filename + ".crpa", iteration, controls_.tracking_definition,
                                 dataset.mf);
                             input_.window.previous_parent_spin_k = tracking_.parent_phi;
                         });
}

void CrpaQsgwExtension::write_checkpoint(const Dataset& dataset, const std::string& filename,
                                         int iteration) const
{
    if (!controls_.track_orbitals) return;
    // QSGW invokes this on its reporting rank; checkpoint serialization has no collectives.
    write_crpa_tracking_checkpoint(filename + ".crpa", iteration, controls_.tracking_definition,
                                   crpa_tracking_energies(dataset.mf), tracking_);
}

void CrpaQsgwExtension::evaluate_final(Dataset& dataset, const LibrpaOptions& options,
                                       int completed_iteration)
{
    emit(dataset, compute_crpa_u_task(dataset, options, input_), completed_iteration,
         completed_iteration, true);
}
}  // namespace librpa_int
