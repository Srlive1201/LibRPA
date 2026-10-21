#include <cmath>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

#include "../mpi/collective_stage.h"
#include "instance_manager.h"
#include "librpa_qsgw.h"
#include "qsgw_helper.h"

struct LibrpaQsgwTask
{
    std::unique_ptr<librpa_int::qsgw::Task> task;
    std::shared_ptr<librpa_int::Dataset> dataset;
};

namespace
{
using namespace librpa_int;
using namespace librpa_int::qsgw;
std::string copied(const char* text, const char* fallback) { return text ? text : fallback; }
Matz decode_matrix(const LibrpaQsgwMatrix& source)
{
    if (source.dimension <= 0 || !source.values_ri ||
        (source.row_major != 0 && source.row_major != 1))
        throw LIBRPA_INVALID_ARGUMENT("QSGW matrix dimensions/storage are invalid");
    Matz matrix(source.dimension, source.dimension, source.row_major ? MAJOR::ROW : MAJOR::COL);
    for (std::size_t i = 0; i < matrix.size(); ++i)
        matrix.ptr()[i] = {source.values_ri[2 * i], source.values_ri[2 * i + 1]};
    return matrix;
}
SpinKMatrixMap decode_matrices(int count, const LibrpaQsgwMatrix* source)
{
    if (count < 0 || (count && !source))
        throw LIBRPA_INVALID_ARGUMENT("QSGW matrix list is incomplete");
    SpinKMatrixMap result;
    for (int i = 0; i < count; ++i)
    {
        const auto& block = source[i];
        if (block.spin < 0 || block.kpoint < 0 ||
            !result[block.spin].emplace(block.kpoint, decode_matrix(block)).second)
            throw LIBRPA_INVALID_ARGUMENT("QSGW matrix repeats or has an invalid state index");
    }
    return result;
}
LocalizedMatrix decode_localized(const LibrpaLocalizedOperator& source, int dimension)
{
    if (source.dimension != dimension || source.n_blocks <= 0 || !source.blocks)
        throw LIBRPA_INVALID_ARGUMENT("QSGW native operator dimensions/blocks are incomplete");
    LocalizedMatrix result;
    for (int i = 0; i < source.n_blocks; ++i)
    {
        const auto& block = source.blocks[i];
        const LibrpaQsgwMatrix view{0, 0, dimension, block.row_major, block.values_ri};
        const Vector3_Order<int> r{block.translation[0], block.translation[1],
                                   block.translation[2]};
        if (!result.emplace(r, decode_matrix(view)).second)
            throw LIBRPA_INVALID_ARGUMENT("QSGW native operator repeats a lattice translation");
    }
    return result;
}
}  // namespace

librpa_int::qsgw::Task& librpa_int::api::qsgw_task_instance(LibrpaQsgwTask* task)
{
    if (!task || !task->task) throw LIBRPA_INVALID_ARGUMENT("Missing QSGW task");
    return *task->task;
}

void librpa_init_qsgw_controls(LibrpaQsgwControls* controls)
{
    if (!controls) throw LIBRPA_INVALID_ARGUMENT("Missing QSGW controls");
    *controls = {1,
                 10,
                 10,
                 0,
                 0.2,
                 2.0 / 27.211386245988,
                 1.0e-4,
                 "none",
                 "reference",
                 "none",
                 "qsgw_checkpoint.bin",
                 0,
                 1,
                 1,
                 0,
                 2147483647,
                 0};
}
LibrpaQsgwTask* librpa_create_qsgw_task(LibrpaHandler* handler, const LibrpaOptions* options,
                                        const LibrpaQsgwControls* controls)
{
    using namespace librpa_int;
    using namespace librpa_int::qsgw;
    auto ds = api::get_dataset_instance(handler);
    TaskControls input;
    collective_all_stage(
        ds->comm_h, "QSGW task controls",
        [&]
        {
            if (!options || !controls)
                throw LIBRPA_INVALID_ARGUMENT("Missing QSGW options/controls");
            input.qsgw_min_iter = controls->min_iterations;
            input.qsgw_max_iter = controls->max_iterations;
            input.qsgw_unoccupied_keep = controls->unoccupied_keep;
            input.qsgw_cut_mode = controls->cut_mode;
            input.qsgw_mixing_beta = controls->mixing_beta;
            input.qsgw_cut_shift_ha = controls->cut_shift_ha;
            input.qsgw_convergence_tolerance_ev = controls->convergence_tolerance_ev;
            input.qsgw_mixer = copied(controls->mixer, "none");
            input.qsgw_headwing_velocity = copied(controls->velocity_mode, "reference");
            input.qsgw_restart_file = copied(controls->restart_file, "none");
            input.qsgw_checkpoint_file = copied(controls->checkpoint_file, "none");
            input.qsgw_write_iteration_matrices = controls->output_matrices != 0;
            input.qsgw_output_qp = controls->output_qp != 0;
            input.qsgw_export_hamiltonian_for_pyatb = controls->export_hamiltonian != 0;
            input.i_state_low = controls->state_low;
            input.i_state_high = controls->state_high;
            input.interband_velocity = controls->interband_velocity != 0;
                if (input.qsgw_min_iter < 1 || input.qsgw_max_iter < input.qsgw_min_iter ||
                input.qsgw_unoccupied_keep < 0 ||
                (input.qsgw_cut_mode != 0 && input.qsgw_cut_mode != 1) ||
                !std::isfinite(input.qsgw_cut_shift_ha) ||
                !std::isfinite(input.qsgw_convergence_tolerance_ev) ||
                input.qsgw_convergence_tolerance_ev <= 0.0 ||
                (input.qsgw_mixer != "none" && input.qsgw_mixer != "linear") ||
                !(input.qsgw_mixing_beta > 0.0 && input.qsgw_mixing_beta <= 1.0) ||
                (input.qsgw_headwing_velocity != "reference" &&
                 input.qsgw_headwing_velocity != "rotate_reference" &&
                 input.qsgw_headwing_velocity != "current_hr"))
                throw LIBRPA_INVALID_ARGUMENT(
                    "Invalid QSGW iteration, mixing or velocity controls");
        });
    auto result = std::make_unique<LibrpaQsgwTask>();
    result->dataset = ds;
    result->task = std::make_unique<Task>(ds, *options, std::move(input));
    return result.release();
}
void librpa_prepare_qsgw_task(LibrpaQsgwTask* task)
{
    librpa_int::api::qsgw_task_instance(task).prepare();
}
void librpa_set_qsgw_input(LibrpaQsgwTask* task, const LibrpaQsgwInput* root_input)
{
    using namespace librpa_int;
    using namespace librpa_int::qsgw;
    auto& engine = api::qsgw_task_instance(task);
    TaskInput input;
    collective_all_stage(
        task->dataset->comm_h, "QSGW numeric input",
        [&]
        {
            if (!task->dataset->comm_h.is_root()) return;
            if (!root_input) throw LIBRPA_INVALID_ARGUMENT("Missing QSGW numeric input on root");
            input.vxc = decode_matrices(root_input->n_vxc, root_input->vxc);
            input.vxc_band = decode_matrices(root_input->n_vxc_band, root_input->vxc_band);
            input.vxc_basis = root_input->vxc_in_ao_basis ? VxcBasis::Nao : VxcBasis::State;
            input.inplane_area_bohr2 = root_input->inplane_area_bohr2;
            input.auxiliary_monopole_norm_squared = root_input->auxiliary_monopole_norm_squared;
            if (const auto* source = root_input->native_operators)
            {
                auto& operators = input.ao_operators.emplace();
                operators.dimension = source->overlap.dimension;
                operators.overlap = decode_localized(source->overlap, operators.dimension);
                for (int axis = 0; axis < 3; ++axis)
                    operators.position[axis] =
                        decode_localized(source->position[axis], operators.dimension);
                if (source->n_spins <= 0 || !source->hamiltonian)
                    throw LIBRPA_INVALID_ARGUMENT("QSGW native Hamiltonians are incomplete");
                for (int spin = 0; spin < source->n_spins; ++spin)
                    operators.native_dft_hamiltonian.push_back(decode_localized(
                        source->hamiltonian[spin],
                        task->dataset->mf.get_n_aos() * task->dataset->mf.get_n_spinor()));
            }
        });
    engine.set_input(std::move(input));
}
void librpa_set_qsgw_output_sink(LibrpaQsgwTask* task, LibrpaQsgwOutputSink sink, void* user_data)
{
    using namespace librpa_int;
    api::qsgw_task_instance(task).set_output(
        [sink, user_data](const qsgw::TaskOutput& event)
        {
            if (!sink) return;
            std::vector<LibrpaLocalizedBlock> blocks;
            LibrpaLocalizedOperator op{};
            if (event.hamiltonian)
            {
                for (const auto& [r, matrix] : *event.hamiltonian)
                    blocks.push_back({{r.x, r.y, r.z},
                                      reinterpret_cast<const double*>(matrix.ptr()),
                                      matrix.major() == MAJOR::ROW});
                op.dimension =
                    event.hamiltonian->empty() ? 0 : event.hamiltonian->begin()->second.nr();
                op.n_blocks = blocks.size();
                op.blocks = blocks.data();
            }
            const LibrpaQsgwOutput output{static_cast<LibrpaQsgwOutputKind>(event.kind),
                                          event.iteration,
                                          event.spin,
                                          event.append,
                                          event.text.c_str(),
                                          event.text.size(),
                                          event.hamiltonian ? &op : nullptr};
            if (sink(&output, user_data) != 0)
                throw LIBRPA_RUNTIME_ERROR("QSGW output sink failed");
        });
}
void librpa_run_qsgw_task(LibrpaQsgwTask* task) { librpa_int::api::qsgw_task_instance(task).run(); }
void librpa_delete_qsgw_task(LibrpaQsgwTask* task) { delete task; }
