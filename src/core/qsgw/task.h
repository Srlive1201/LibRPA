#pragma once
#include <functional>
#include <memory>
#include <optional>
#include <string>

#include "task_context.h"
#include "../meanfield_mpi.h"
#include "headwing_velocity.h"
#include "librpa_options.h"
#include "task_extension.h"
#include "vxc_projection.h"

namespace librpa_int::qsgw
{
struct TaskControls
{
    std::string qsgw_mixer = "none";
    double qsgw_mixing_beta = 0.2;
    int qsgw_min_iter = 1, qsgw_max_iter = 10;
    int qsgw_unoccupied_keep = 10, qsgw_cut_mode = 0;
    double qsgw_cut_shift_ha = 2.0 / 27.211386245988;
    bool qsgw_write_iteration_matrices = false, qsgw_output_qp = true;
    bool qsgw_export_hamiltonian_for_pyatb = true;
    std::string qsgw_headwing_velocity = "reference";
    std::string qsgw_restart_file = "none", qsgw_checkpoint_file = "qsgw_checkpoint.bin";
    double qsgw_convergence_tolerance_ev = 1.0e-4;
    int i_state_low = 0, i_state_high = 2147483647;
    bool interband_velocity = false, expand_ibz = false;
};
enum class OutputKind
{
    Iterations,
    Eigenvalues,
    Matrices,
    Quasiparticles,
    EnergyQp,
    GwBand,
    KsBand,
    ExxBand,
    QsgwBand,
    Hamiltonian,
    HeadHamiltonian
};
struct TaskOutput
{
    OutputKind kind;
    int iteration = 0, spin = 0;
    bool append = false;
    std::string text;
    const LocalizedMatrix* hamiltonian = nullptr;
};
struct TaskInput
{
    SpinKMatrixMap vxc, vxc_band;
    VxcBasis vxc_basis = VxcBasis::State;
    std::optional<HeadwingAoOperators> ao_operators;
    double inplane_area_bohr2 = 0.0, auxiliary_monopole_norm_squared = 0.0;
};
class Task
{
public:
    Task(std::shared_ptr<TaskContext> dataset, const LibrpaOptions&, TaskControls);
    // Validate controls before ordinary numerical input is installed.
    void prepare();
    void set_input(TaskInput input);
    void set_output(std::function<void(const TaskOutput&)> sink);
    void set_extension(std::shared_ptr<TaskExtension> extension);
    void run();

private:
    void output_text(OutputKind, const std::string&, int iteration = 0, int spin = 0,
                     bool append = false) const;
    void output_hamiltonian(OutputKind, const LocalizedMatrix&, int iteration, int spin) const;
    std::shared_ptr<TaskContext> dataset_;
    LibrpaOptions options_;
    TaskControls controls_;
    TaskInput input_;
    std::function<void(const TaskOutput&)> sink_;
    std::shared_ptr<TaskExtension> extension_;
    std::optional<MeanField> ibz_reference_;
    std::optional<SymmetryExpandedMeanField> expanded_reference_;
    bool prepared_ = false, input_set_ = false, ran_ = false;
};
}  // namespace librpa_int::qsgw
