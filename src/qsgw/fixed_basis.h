#pragma once

#include "matrix_map.h"

#include "../core/meanfield.h"
#include "../math/complexmatrix.h"
#include "../mpi/base_mpi.h"

#include <map>
#include <vector>

namespace librpa_int
{
namespace qsgw
{

using VelocityMatrix =
    std::vector<std::vector<std::vector<ComplexMatrix>>>;

using MeanFieldEigenvectorMap =
    std::map<int, std::map<int, std::map<int, ComplexMatrix>>>;

// Temporarily expose the immutable reference WFC through a live MeanField so
// the unmodified upstream k-grid projection routines can be reused. All other
// live mean-field state remains untouched and the live WFC is restored on exit.
class ScopedReferenceEigenvectors
{
public:
    ScopedReferenceEigenvectors(MeanField& live,
                                const MeanField& reference);
    ~ScopedReferenceEigenvectors() noexcept;

    ScopedReferenceEigenvectors(const ScopedReferenceEigenvectors&) = delete;
    ScopedReferenceEigenvectors& operator=(
        const ScopedReferenceEigenvectors&) = delete;
    ScopedReferenceEigenvectors(ScopedReferenceEigenvectors&&) = delete;
    ScopedReferenceEigenvectors& operator=(
        ScopedReferenceEigenvectors&&) = delete;

private:
    MeanField& live_;
    MeanFieldEigenvectorMap live_eigenvectors_;
};

struct FixedBasisDiagonalizationResult
{
    SpinKMatrixMap unitary;
};

struct VelocityBasisAlignmentResult
{
    double maximum_relative_wfc_residual = 0.0;
    double maximum_unitarity_residual = 0.0;
    double maximum_raw_relative_wfc_residual = 0.0;
    double maximum_raw_unitarity_residual = 0.0;
    double maximum_unitary_projection_correction = 0.0;
    double maximum_basis_inverse_residual = 0.0;
    double maximum_basis_condition_estimate = 0.0;
    double maximum_transform_deviation_from_identity = 0.0;
};

void prepare_fhi_aims_interband_velocity(
    VelocityMatrix& velocity,
    const MeanField& reference);

VelocityBasisAlignmentResult align_velocity_to_reference_wfc(
    const MeanField& velocity_basis,
    const MeanField& reference,
    VelocityMatrix& velocity,
    double relative_tolerance = 1.0e-8);

VelocityBasisAlignmentResult align_distributed_velocity_to_reference_wfc(
    const MeanField& velocity_basis,
    const MeanField& reference,
    VelocityMatrix& velocity,
    const MpiCommHandler& communicator,
    double relative_tolerance = 1.0e-8);

FixedBasisDiagonalizationResult diagonalize_in_reference_basis(
    MeanField& live,
    const MeanField& reference,
    const SpinKMatrixMap& hamiltonian);

} // namespace qsgw
} // namespace librpa_int
