#pragma once

#include "../../src/core/qsgw/hamiltonian_mixing.h"
#include "../../src/core/qsgw/headwing_velocity.h"
#include "../../src/mpi/base_mpi.h"
#include "iteration_trace.h"

namespace librpa_int
{
namespace qsgw
{

// The checkpoint remains in the immutable input KS basis. Native producer
// inputs are still required to rebuild RI, Coulomb and velocity objects.
// With current_hr, ao_operators includes the immutable native DFT H(R) baseline
// (in Ha) used by both the exported Hamiltonian and its velocity derivatives.
struct CheckpointContext
{
    const MeanField& reference;
    const MeanField* band_reference;
    const std::vector<Vector3_Order<double>>& kpoints;
    const std::vector<Vector3_Order<double>>& band_kpoints;
    const std::vector<double>& k_weights;
    const SpinKMatrixMap& dft_vxc;
    const SpinKMatrixMap& dft_vxc_band;
    int cut_mode;
    int unoccupied_keep;
    double cut_shift_ha;
    double mixing_beta;
    int nfreq;
    int head_mode;
    std::string velocity_mode = "reference";
    const HeadwingAoOperators* ao_operators = nullptr;
    const VelocityMatrix* reference_velocity = nullptr;
    const Vector3_Order<int>* velocity_period = nullptr;
    const Matrix3* velocity_lattice = nullptr;
};

// Root-only file operations; callers propagate failures collectively.
void write_qsgw_checkpoint(const std::string& filename, const CheckpointContext& context,
                           const IterationSummary& summary, const MeanField& live,
                           const MeanField* band_live, const SpinKMatrixMap& hamiltonian,
                           const SpinKMatrixMap& band_hamiltonian);

IterationSummary read_qsgw_checkpoint(const std::string& filename, const CheckpointContext& context,
                                      MeanField& live, MeanField* band_live,
                                      SpinKMatrixMap& hamiltonian,
                                      SpinKMatrixMap& band_hamiltonian);

// QSGW currently requires complete replicated wavefunctions on every rank.
void broadcast_checkpoint_meanfield(MeanField& live, const MpiCommHandler& communicator);

}  // namespace qsgw
}  // namespace librpa_int
