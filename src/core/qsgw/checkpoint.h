#pragma once
#include <string>

#include "distributed_matrix.h"
#include "headwing_velocity.h"
#include "iteration_trace.h"
namespace librpa_int::qsgw
{
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

IterationSummary read_qsgw_checkpoint(const std::string&, const CheckpointContext&, MeanField&,
                                      MeanField*, SpinKMatrixMap&, SpinKMatrixMap&);
void write_qsgw_checkpoint(const std::string&, const CheckpointContext&, const IterationSummary&,
                           const MeanField&, const MeanField*, const SpinKMatrixMap&,
                           const SpinKMatrixMap&);
void broadcast_checkpoint_meanfield(MeanField&, const MpiCommHandler&);
}  // namespace librpa_int::qsgw
