#include "fixed_basis.h"

#include "../math/utils_matrix_mpi.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace librpa_int
{
namespace qsgw
{
namespace
{

constexpr double hermitian_tolerance = 1.0e-10;

void require_same_meanfield_shape(const MeanField& live,
                                  const MeanField& reference)
{
    if (live.get_n_spins() != reference.get_n_spins() ||
        live.get_n_kpoints() != reference.get_n_kpoints() ||
        live.get_n_bands() != reference.get_n_bands() ||
        live.get_n_aos() != reference.get_n_aos() ||
        live.get_n_spinor() != reference.get_n_spinor())
    {
        throw std::invalid_argument(
            "QSGW live and reference mean fields have different shapes");
    }
}

void require_same_wfc_distribution(const MeanField& live,
                                   const MeanField& reference)
{
    const auto& live_wfc = live.get_eigenvectors();
    const auto& reference_wfc = reference.get_eigenvectors();
    if (live_wfc.size() != reference_wfc.size())
    {
        throw std::invalid_argument(
            "QSGW live and reference WFC distributions differ");
    }
    for (const auto& [spin, reference_spin] : reference_wfc)
    {
        const auto live_spin_it = live_wfc.find(spin);
        if (live_spin_it == live_wfc.end() ||
            live_spin_it->second.size() != reference_spin.size())
        {
            throw std::invalid_argument(
                "QSGW live and reference WFC distributions differ");
        }
        for (const auto& [spinor, reference_spinor] : reference_spin)
        {
            const auto live_spinor_it = live_spin_it->second.find(spinor);
            if (live_spinor_it == live_spin_it->second.end() ||
                live_spinor_it->second.size() != reference_spinor.size())
            {
                throw std::invalid_argument(
                    "QSGW live and reference WFC distributions differ");
            }
            for (const auto& [kpoint, reference_block] : reference_spinor)
            {
                const auto live_block_it =
                    live_spinor_it->second.find(kpoint);
                if (live_block_it == live_spinor_it->second.end() ||
                    live_block_it->second.nr != reference_block.nr ||
                    live_block_it->second.nc != reference_block.nc)
                {
                    throw std::invalid_argument(
                        "QSGW live and reference WFC distributions differ");
                }
            }
        }
    }
}

void require_valid_wfc(const MeanField& meanfield, const char* label)
{
    for (int spin = 0; spin < meanfield.get_n_spins(); ++spin)
    {
        for (int spinor = 0; spinor < meanfield.get_n_spinor(); ++spinor)
        {
            for (int kpoint = 0; kpoint < meanfield.get_n_kpoints(); ++kpoint)
            {
                const ComplexMatrix* block =
                    meanfield.find_wfc(spin, spinor, kpoint);
                if (block == nullptr ||
                    block->nr != meanfield.get_n_bands() ||
                    block->nc != meanfield.get_n_aos())
                {
                    throw std::invalid_argument(
                        std::string("QSGW ") + label +
                        " wavefunction map is incomplete or has an invalid shape");
                }
                for (int index = 0; index < block->size; ++index)
                {
                    if (!std::isfinite(block->c[index].real()) ||
                        !std::isfinite(block->c[index].imag()))
                    {
                        throw std::invalid_argument(
                            std::string("QSGW ") + label +
                            " wavefunction contains non-finite data");
                    }
                }
            }
        }
    }
}

void require_valid_hamiltonian(const SpinKMatrixMap& hamiltonian,
                               const MeanField& reference)
{
    if (static_cast<int>(hamiltonian.size()) != reference.get_n_spins())
    {
        throw std::invalid_argument("QSGW Hamiltonian spin map is incomplete");
    }
    const int dimension = reference.get_n_bands();
    for (int spin = 0; spin < reference.get_n_spins(); ++spin)
    {
        const auto spin_it = hamiltonian.find(spin);
        if (spin_it == hamiltonian.end() ||
            static_cast<int>(spin_it->second.size()) !=
                reference.get_n_kpoints())
        {
            throw std::invalid_argument("QSGW Hamiltonian k-point map is incomplete");
        }
        for (int kpoint = 0; kpoint < reference.get_n_kpoints(); ++kpoint)
        {
            const auto kpoint_it = spin_it->second.find(kpoint);
            if (kpoint_it == spin_it->second.end())
            {
                throw std::invalid_argument(
                    "QSGW Hamiltonian k-point map is incomplete");
            }
            const Matz& matrix = kpoint_it->second;
            if (matrix.nr() != dimension || matrix.nc() != dimension)
            {
                throw std::invalid_argument(
                    "QSGW Hamiltonian matrix has an invalid shape");
            }
            for (int row = 0; row < dimension; ++row)
            {
                for (int column = 0; column < dimension; ++column)
                {
                    const cplxdb value = matrix(row, column);
                    if (!std::isfinite(value.real()) ||
                        !std::isfinite(value.imag()))
                    {
                        throw std::invalid_argument(
                            "QSGW Hamiltonian contains non-finite data");
                    }
                    if (std::abs(value - std::conj(matrix(column, row))) >
                        hermitian_tolerance)
                    {
                        throw std::invalid_argument(
                            "QSGW Hamiltonian is not Hermitian");
                    }
                }
            }
        }
    }
}

void require_valid_velocity(const VelocityMatrix& velocity,
                            const MeanField& reference,
                            const char* label)
{
    if (static_cast<int>(velocity.size()) != reference.get_n_spins())
    {
        throw std::invalid_argument(
            std::string("QSGW ") + label + " velocity spin map is incomplete");
    }
    const int dimension = reference.get_n_bands();
    for (int spin = 0; spin < reference.get_n_spins(); ++spin)
    {
        if (static_cast<int>(velocity[spin].size()) !=
            reference.get_n_kpoints())
        {
            throw std::invalid_argument(
                std::string("QSGW ") + label +
                " velocity k-point map is incomplete");
        }
        for (int kpoint = 0; kpoint < reference.get_n_kpoints(); ++kpoint)
        {
            if (velocity[spin][kpoint].size() != 3)
            {
                throw std::invalid_argument(
                    std::string("QSGW ") + label +
                    " velocity must contain three Cartesian components");
            }
            for (const ComplexMatrix& component : velocity[spin][kpoint])
            {
                if (component.nr != dimension || component.nc != dimension)
                {
                    throw std::invalid_argument(
                        std::string("QSGW ") + label +
                        " velocity matrix has an invalid shape");
                }
                for (int index = 0; index < component.size; ++index)
                {
                    if (!std::isfinite(component.c[index].real()) ||
                        !std::isfinite(component.c[index].imag()))
                    {
                        throw std::invalid_argument(
                            std::string("QSGW ") + label +
                            " velocity contains non-finite data");
                    }
                }
            }
        }
    }
}

Matz collect_reference_wfc_rows(const MeanField& reference,
                                const int spin,
                                const int kpoint,
                                const MAJOR major)
{
    const int dimension = reference.get_n_bands();
    const int spinor_count = reference.get_n_spinor();
    const int ao_count = reference.get_n_aos();
    Matz rows(dimension, ao_count * spinor_count, major);
    for (int band = 0; band < dimension; ++band)
    {
        for (int spinor = 0; spinor < spinor_count; ++spinor)
        {
            const ComplexMatrix& block =
                reference.get_eigenvectors().at(spin).at(spinor).at(kpoint);
            for (int ao = 0; ao < ao_count; ++ao)
            {
                rows(band, ao * spinor_count + spinor) = block(band, ao);
            }
        }
    }
    return rows;
}

void store_wfc_rows(MeanField& target,
                    const int spin,
                    const int kpoint,
                    const Matz& rows)
{
    for (int band = 0; band < target.get_n_bands(); ++band)
    {
        for (int spinor = 0; spinor < target.get_n_spinor(); ++spinor)
        {
            ComplexMatrix& block =
                target.get_eigenvectors().at(spin).at(spinor).at(kpoint);
            for (int ao = 0; ao < target.get_n_aos(); ++ao)
            {
                block(band, ao) =
                    rows(band, ao * target.get_n_spinor() + spinor);
            }
        }
    }
}

ComplexMatrix to_complex_matrix(const Matz& input)
{
    ComplexMatrix output(input.nr(), input.nc());
    for (int row = 0; row < input.nr(); ++row)
    {
        for (int column = 0; column < input.nc(); ++column)
        {
            output(row, column) = input(row, column);
        }
    }
    return output;
}

double frobenius_norm(const Matz& matrix)
{
    double norm_squared = 0.0;
    for (int row = 0; row < matrix.nr(); ++row)
    {
        for (int column = 0; column < matrix.nc(); ++column)
        {
            norm_squared += std::norm(matrix(row, column));
        }
    }
    return std::sqrt(norm_squared);
}

void require_finite_matrix(const Matz& matrix, const char* label)
{
    for (int row = 0; row < matrix.nr(); ++row)
    {
        for (int column = 0; column < matrix.nc(); ++column)
        {
            const cplxdb value = matrix(row, column);
            if (!std::isfinite(value.real()) ||
                !std::isfinite(value.imag()))
            {
                throw std::invalid_argument(
                    std::string("QSGW ") + label +
                    " contains non-finite data");
            }
        }
    }
}

Matz invert_complete_wfc_basis(
    const Matz& matrix,
    const double relative_tolerance,
    VelocityBasisAlignmentResult& diagnostics)
{
    if (matrix.nr() != matrix.nc())
    {
        throw std::invalid_argument(
            "QSGW velocity alignment requires a complete square WFC basis");
    }
    const int dimension = matrix.nr();
    Matz result = matrix.copy();
    std::vector<int> pivots(static_cast<std::size_t>(dimension));
    std::vector<cplxdb> work(
        static_cast<std::size_t>(std::max(1, dimension)));
    int info = 0;
    if (result.is_row_major())
    {
        LapackConnector::getrf(
            dimension, dimension, result.ptr(), dimension,
            pivots.data(), info);
    }
    else
    {
        LapackConnector::getrf_f(
            dimension, dimension, result.ptr(), dimension,
            pivots.data(), info);
    }
    if (info != 0)
    {
        throw std::invalid_argument(
            "QSGW velocity reference WFC basis is singular");
    }
    if (result.is_row_major())
    {
        LapackConnector::getri(
            dimension, result.ptr(), dimension, pivots.data(),
            work.data(), static_cast<int>(work.size()), info);
    }
    else
    {
        LapackConnector::getri_f(
            dimension, result.ptr(), dimension, pivots.data(),
            work.data(), static_cast<int>(work.size()), info);
    }
    if (info != 0)
    {
        throw std::invalid_argument(
            "QSGW velocity reference WFC basis inversion failed");
    }
    require_finite_matrix(result, "velocity reference WFC inverse");

    const Matz product = matrix * result;
    double residual_squared = 0.0;
    for (int row = 0; row < dimension; ++row)
    {
        for (int column = 0; column < dimension; ++column)
        {
            const cplxdb expected = row == column ? 1.0 : 0.0;
            residual_squared +=
                std::norm(product(row, column) - expected);
        }
    }
    const double inverse_residual = std::sqrt(
        residual_squared / std::max(1, dimension));
    const double condition_estimate =
        frobenius_norm(matrix) * frobenius_norm(result);
    diagnostics.maximum_basis_inverse_residual = std::max(
        diagnostics.maximum_basis_inverse_residual,
        inverse_residual);
    diagnostics.maximum_basis_condition_estimate = std::max(
        diagnostics.maximum_basis_condition_estimate,
        condition_estimate);
    if (!std::isfinite(inverse_residual) ||
        inverse_residual > relative_tolerance ||
        !std::isfinite(condition_estimate) ||
        condition_estimate > 1.0 / relative_tolerance)
    {
        throw std::invalid_argument(
            "QSGW velocity reference WFC basis is ill-conditioned");
    }
    return result;
}

void diagonalize_hermitian(const Matz& matrix,
                           std::vector<double>& eigenvalues,
                           Matz& unitary)
{
    const int dimension = matrix.nr();
    const int block_size = LapackConnector::ilaenv(
        1, "zheev", "VU", dimension, -1, -1, -1);
    const int work_size = std::max(1, dimension * (block_size + 1));
    const int real_work_size = std::max(1, 3 * dimension - 2);

    unitary = matrix.copy();
    eigenvalues.assign(dimension, 0.0);
    std::vector<cplxdb> work(work_size);
    std::vector<double> real_work(real_work_size);
    int info = 0;
    if (unitary.is_row_major())
    {
        LapackConnector::heev(
            'V', 'U', dimension, unitary.ptr(), dimension,
            eigenvalues.data(), work.data(), work_size,
            real_work.data(), info);
    }
    else
    {
        LapackConnector::heev_f(
            'V', 'U', dimension, unitary.ptr(), dimension,
            eigenvalues.data(), work.data(), work_size,
            real_work.data(), info);
    }
    if (info != 0)
    {
        throw std::runtime_error(
            "QSGW fixed-basis Hamiltonian diagonalization failed");
    }
}

double relative_reconstruction_residual(const Matz& source,
                                        const Matz& transform,
                                        const Matz& reference)
{
    const Matz reconstructed = transform * reference;
    double residual_squared = 0.0;
    double source_norm_squared = 0.0;
    double reconstructed_norm_squared = 0.0;
    for (int row = 0; row < source.nr(); ++row)
    {
        for (int column = 0; column < source.nc(); ++column)
        {
            residual_squared += std::norm(
                source(row, column) - reconstructed(row, column));
            source_norm_squared += std::norm(source(row, column));
            reconstructed_norm_squared +=
                std::norm(reconstructed(row, column));
        }
    }
    return std::sqrt(
        residual_squared /
        std::max(source_norm_squared, reconstructed_norm_squared));
}

double maximum_unitarity_residual(const Matz& transform)
{
    const Matz product = transform * transpose(transform, true);
    double result = 0.0;
    for (int row = 0; row < transform.nr(); ++row)
    {
        for (int column = 0; column < transform.nc(); ++column)
        {
            const cplxdb expected = row == column ? 1.0 : 0.0;
            result = std::max(
                result, std::abs(product(row, column) - expected));
        }
    }
    return result;
}

Matz project_to_nearest_unitary(const Matz& transform)
{
    const int dimension = transform.nr();
    const Matz gram = transpose(transform, true) * transform;
    std::vector<double> eigenvalues;
    Matz eigenvectors;
    diagonalize_hermitian(gram, eigenvalues, eigenvectors);
    for (const double eigenvalue : eigenvalues)
    {
        if (!(eigenvalue > 0.0) || !std::isfinite(eigenvalue))
        {
            throw std::invalid_argument(
                "QSGW velocity WFC basis transform has a non-positive singular value");
        }
    }

    Matz inverse_sqrt(dimension, dimension, transform.major());
    for (int row = 0; row < dimension; ++row)
    {
        for (int column = 0; column < dimension; ++column)
        {
            cplxdb value = 0.0;
            for (int state = 0; state < dimension; ++state)
            {
                value += eigenvectors(row, state) *
                         (1.0 / std::sqrt(eigenvalues[state])) *
                         std::conj(eigenvectors(column, state));
            }
            inverse_sqrt(row, column) = value;
        }
    }
    Matz result = transform * inverse_sqrt;
    require_finite_matrix(result, "projected velocity WFC basis transform");
    return result;
}

} // namespace

ScopedReferenceEigenvectors::ScopedReferenceEigenvectors(
    MeanField& live,
    const MeanField& reference)
    : live_(live),
      live_eigenvectors_(reference.get_eigenvectors())
{
    if (&live == &reference)
    {
        throw std::invalid_argument(
            "QSGW live and reference mean fields must be distinct");
    }
    require_same_meanfield_shape(live, reference);
    require_same_wfc_distribution(live, reference);
    live_.get_eigenvectors().swap(live_eigenvectors_);
}

ScopedReferenceEigenvectors::~ScopedReferenceEigenvectors() noexcept
{
    live_.get_eigenvectors().swap(live_eigenvectors_);
}

void prepare_fhi_aims_interband_velocity(
    VelocityMatrix& velocity,
    const MeanField& reference)
{
    require_valid_velocity(velocity, reference, "FHI-aims source");
    VelocityMatrix prepared = velocity;
    for (auto& spin : prepared)
    {
        for (auto& kpoint : spin)
        {
            for (ComplexMatrix& component : kpoint)
            {
                for (int row = 0; row < component.nr; ++row)
                {
                    for (int column = row + 1;
                         column < component.nc; ++column)
                    {
                        if (std::abs(component(row, column) -
                                     std::conj(component(column, row))) >
                            hermitian_tolerance)
                        {
                            throw std::invalid_argument(
                                "QSGW FHI-aims interband velocity is not Hermitian");
                        }
                    }
                    component(row, row) = 0.0;
                }
            }
        }
    }
    velocity = std::move(prepared);
}

VelocityBasisAlignmentResult align_velocity_to_reference_wfc(
    const MeanField& velocity_basis,
    const MeanField& reference,
    VelocityMatrix& velocity,
    const double relative_tolerance)
{
    if (!(relative_tolerance > 0.0) ||
        !std::isfinite(relative_tolerance))
    {
        throw std::invalid_argument(
            "QSGW velocity-basis tolerance must be finite and positive");
    }
    require_same_meanfield_shape(velocity_basis, reference);
    require_valid_wfc(velocity_basis, "velocity-basis");
    require_valid_wfc(reference, "reference");
    require_valid_velocity(velocity, reference, "source");

    VelocityMatrix aligned = velocity;
    VelocityBasisAlignmentResult result;
    const int dimension = reference.get_n_bands();
    const int coefficient_dimension =
        reference.get_n_aos() * reference.get_n_spinor();
    if (dimension != coefficient_dimension)
    {
        throw std::invalid_argument(
            "QSGW velocity alignment requires a complete square WFC basis");
    }
    for (int spin = 0; spin < reference.get_n_spins(); ++spin)
    {
        for (int kpoint = 0; kpoint < reference.get_n_kpoints(); ++kpoint)
        {
            const Matz source_rows = collect_reference_wfc_rows(
                velocity_basis, spin, kpoint, MAJOR::ROW);
            const Matz reference_rows = collect_reference_wfc_rows(
                reference, spin, kpoint, MAJOR::ROW);
            const Matz reference_inverse = invert_complete_wfc_basis(
                reference_rows, relative_tolerance, result);
            const Matz raw_transform = source_rows * reference_inverse;
            require_finite_matrix(
                raw_transform, "raw velocity WFC basis transform");

            const double raw_relative_residual =
                relative_reconstruction_residual(
                    source_rows, raw_transform, reference_rows);
            result.maximum_raw_relative_wfc_residual = std::max(
                result.maximum_raw_relative_wfc_residual,
                raw_relative_residual);
            if (!std::isfinite(raw_relative_residual) ||
                raw_relative_residual > relative_tolerance)
            {
                throw std::invalid_argument(
                    "QSGW raw velocity WFC basis transform does not reconstruct the source basis");
            }

            const double raw_unitarity_residual =
                maximum_unitarity_residual(raw_transform);
            result.maximum_raw_unitarity_residual = std::max(
                result.maximum_raw_unitarity_residual,
                raw_unitarity_residual);
            if (!std::isfinite(raw_unitarity_residual) ||
                raw_unitarity_residual > relative_tolerance)
            {
                throw std::invalid_argument(
                    "QSGW raw velocity WFC basis transform is not unitary within input precision");
            }

            const Matz transform =
                project_to_nearest_unitary(raw_transform);
            const double relative_residual =
                relative_reconstruction_residual(
                    source_rows, transform, reference_rows);
            result.maximum_relative_wfc_residual = std::max(
                result.maximum_relative_wfc_residual,
                relative_residual);
            if (!std::isfinite(relative_residual) ||
                relative_residual > relative_tolerance)
            {
                throw std::invalid_argument(
                    "QSGW projected velocity WFC basis transform does not reconstruct the source basis");
            }

            const double unitary_residual =
                maximum_unitarity_residual(transform);
            result.maximum_unitarity_residual = std::max(
                result.maximum_unitarity_residual,
                unitary_residual);
            if (!std::isfinite(unitary_residual) ||
                unitary_residual > relative_tolerance)
            {
                throw std::invalid_argument(
                    "QSGW projected velocity WFC basis transform is not unitary");
            }

            for (int row = 0; row < dimension; ++row)
            {
                for (int column = 0; column < dimension; ++column)
                {
                    const cplxdb expected = row == column ? 1.0 : 0.0;
                    result.maximum_unitary_projection_correction =
                        std::max(
                            result.maximum_unitary_projection_correction,
                            std::abs(
                                transform(row, column) -
                                raw_transform(row, column)));
                    result.maximum_transform_deviation_from_identity =
                        std::max(
                            result.maximum_transform_deviation_from_identity,
                            std::abs(transform(row, column) - expected));
                }
            }
            const ComplexMatrix transform_matrix =
                to_complex_matrix(transform);
            const ComplexMatrix transform_transpose =
                transpose(transform_matrix, false);
            const ComplexMatrix transform_conjugate =
                librpa_int::conj(transform_matrix);
            for (int direction = 0; direction < 3; ++direction)
            {
                aligned[spin][kpoint][direction] =
                    transform_transpose *
                    velocity[spin][kpoint][direction] *
                    transform_conjugate;
            }
        }
    }
    velocity = std::move(aligned);
    return result;
}

VelocityBasisAlignmentResult align_distributed_velocity_to_reference_wfc(
    const MeanField& velocity_basis,
    const MeanField& reference,
    VelocityMatrix& velocity,
    const MpiCommHandler& communicator,
    const double relative_tolerance)
{
    if (!communicator.is_initialized())
    {
        throw std::invalid_argument(
            "QSGW velocity-basis communicator is not initialized");
    }
    if (!(relative_tolerance > 0.0) ||
        !std::isfinite(relative_tolerance))
    {
        throw std::invalid_argument(
            "QSGW velocity-basis tolerance must be finite and positive");
    }
    require_same_meanfield_shape(velocity_basis, reference);

    MeanField replicated_basis = velocity_basis;
    const int expected_rows = reference.get_n_bands();
    const int expected_columns = reference.get_n_aos();
    for (int spin = 0; spin < reference.get_n_spins(); ++spin)
    {
        for (int spinor = 0;
             spinor < reference.get_n_spinor(); ++spinor)
        {
            for (int kpoint = 0;
                 kpoint < reference.get_n_kpoints(); ++kpoint)
            {
                const ComplexMatrix* local =
                    velocity_basis.find_wfc(spin, spinor, kpoint);
                const int owns_local = local == nullptr ? 0 : 1;
                const int bad_shape_local =
                    local != nullptr &&
                    (local->nr != expected_rows ||
                     local->nc != expected_columns)
                        ? 1
                        : 0;
                int bad_shape_global = 0;
                communicator.allreduce(
                    &bad_shape_local, &bad_shape_global, 1, MPI_MAX);
                if (bad_shape_global != 0)
                {
                    throw std::invalid_argument(
                        "QSGW distributed velocity-basis wavefunction has an invalid shape");
                }

                const int owner = find_mpi_owner_rank(
                    owns_local != 0, communicator.comm);
                if (owner < 0)
                {
                    throw std::invalid_argument(
                        "QSGW distributed velocity-basis wavefunction is missing on every MPI rank");
                }

                ComplexMatrix canonical;
                if (communicator.myid == owner)
                {
                    canonical = *local;
                }
                broadcast_ComplexMatrix(
                    canonical, owner, communicator.comm);

                int nonfinite_local = 0;
                double difference_local = 0.0;
                double scale_local = 0.0;
                if (local != nullptr)
                {
                    for (int index = 0; index < local->size; ++index)
                    {
                        const cplxdb value = local->c[index];
                        const cplxdb reference_value = canonical.c[index];
                        if (!std::isfinite(value.real()) ||
                            !std::isfinite(value.imag()) ||
                            !std::isfinite(reference_value.real()) ||
                            !std::isfinite(reference_value.imag()))
                        {
                            nonfinite_local = 1;
                            continue;
                        }
                        difference_local = std::max(
                            difference_local,
                            std::abs(value - reference_value));
                        scale_local = std::max(
                            scale_local, std::abs(reference_value));
                    }
                }
                int nonfinite_global = 0;
                double difference_global = 0.0;
                double scale_global = 0.0;
                communicator.allreduce(
                    &nonfinite_local, &nonfinite_global, 1, MPI_MAX);
                communicator.allreduce(
                    &difference_local, &difference_global, 1, MPI_MAX);
                communicator.allreduce(
                    &scale_local, &scale_global, 1, MPI_MAX);
                if (nonfinite_global != 0 ||
                    difference_global >
                        relative_tolerance * std::max(1.0, scale_global))
                {
                    throw std::invalid_argument(
                        "QSGW distributed velocity-basis wavefunction copies disagree across MPI ranks");
                }

                replicated_basis.get_eigenvectors()[spin][spinor][kpoint] =
                    std::move(canonical);
            }
        }
    }

    return align_velocity_to_reference_wfc(
        replicated_basis, reference, velocity, relative_tolerance);
}

FixedBasisDiagonalizationResult diagonalize_in_reference_basis(
    MeanField& live,
    const MeanField& reference,
    const SpinKMatrixMap& hamiltonian)
{
    if (&live == &reference)
    {
        throw std::invalid_argument(
            "QSGW live and reference mean fields must be distinct objects");
    }
    require_same_meanfield_shape(live, reference);
    require_valid_wfc(reference, "reference");
    require_valid_wfc(live, "live");
    require_valid_hamiltonian(hamiltonian, reference);

    MeanField next_live = live;

    FixedBasisDiagonalizationResult result;
    for (int spin = 0; spin < reference.get_n_spins(); ++spin)
    {
        for (int kpoint = 0; kpoint < reference.get_n_kpoints(); ++kpoint)
        {
            const Matz& matrix = hamiltonian.at(spin).at(kpoint);
            std::vector<double> eigenvalues;
            Matz unitary;
            diagonalize_hermitian(matrix, eigenvalues, unitary);
            result.unitary[spin][kpoint] = unitary;

            for (int band = 0; band < reference.get_n_bands(); ++band)
            {
                next_live.get_eigenvals()[spin](kpoint, band) =
                    eigenvalues.at(band);
            }

            const Matz reference_rows = collect_reference_wfc_rows(
                reference, spin, kpoint, unitary.major());
            const Matz rotated_rows =
                transpose(unitary, false) * reference_rows;
            store_wfc_rows(next_live, spin, kpoint, rotated_rows);
        }
    }

    live = std::move(next_live);
    return result;
}

} // namespace qsgw
} // namespace librpa_int
