#ifdef NDEBUG
#undef NDEBUG
#endif

#include "../qsgw/fixed_basis.h"

#include <cassert>
#include <cmath>
#include <complex>
#include <stdexcept>
#include <utility>
#include <vector>

#include "../io/global_io.h"
#include "../mpi/global_mpi.h"

using librpa_int::ComplexMatrix;
using librpa_int::MeanField;
using librpa_int::cplxdb;
using librpa_int::global::finalize_global_io;
using librpa_int::global::finalize_global_mpi;
using librpa_int::global::init_global_io;
using librpa_int::global::init_global_mpi;
using librpa_int::global::mpi_comm_global_h;
using librpa_int::global::myid_global;
using librpa_int::global::size_global;
using librpa_int::qsgw::VelocityMatrix;
using librpa_int::qsgw::align_distributed_velocity_to_reference_wfc;

namespace
{

void initialize_velocity(VelocityMatrix& velocity,
                         const int kpoints,
                         const int states)
{
    velocity.assign(1, {});
    velocity[0].assign(kpoints, {});
    for (int kpoint = 0; kpoint < kpoints; ++kpoint)
    {
        velocity[0][kpoint].assign(
            3, ComplexMatrix(states, states));
    }
}

void assert_close(const cplxdb actual,
                  const cplxdb expected,
                  const double tolerance = 1.0e-12)
{
    assert(std::abs(actual - expected) < tolerance);
}

void assert_velocity_equal(const VelocityMatrix& actual,
                           const VelocityMatrix& expected)
{
    assert(actual.size() == expected.size());
    for (std::size_t spin = 0; spin < actual.size(); ++spin)
    {
        assert(actual[spin].size() == expected[spin].size());
        for (std::size_t kpoint = 0;
             kpoint < actual[spin].size(); ++kpoint)
        {
            assert(actual[spin][kpoint].size() ==
                   expected[spin][kpoint].size());
            for (std::size_t direction = 0;
                 direction < actual[spin][kpoint].size(); ++direction)
            {
                const auto& lhs = actual[spin][kpoint][direction];
                const auto& rhs = expected[spin][kpoint][direction];
                assert(lhs.nr == rhs.nr);
                assert(lhs.nc == rhs.nc);
                for (int row = 0; row < lhs.nr; ++row)
                {
                    for (int column = 0; column < lhs.nc; ++column)
                    {
                        assert_close(lhs(row, column),
                                     rhs(row, column));
                    }
                }
            }
        }
    }
}

void test_distributed_velocity_basis_alignment()
{
    assert(size_global == 4);
    constexpr int states = 2;
    const int kpoints = size_global;
    MeanField reference(1, kpoints, states, states, 1);
    MeanField distributed_basis(1, kpoints, states, states, 1);
    VelocityMatrix expected_velocity;
    VelocityMatrix source_velocity;
    initialize_velocity(expected_velocity, kpoints, states);

    for (int kpoint = 0; kpoint < kpoints; ++kpoint)
    {
        ComplexMatrix identity(states, states);
        identity.set_as_identity_matrix();
        reference.get_eigenvectors()[0][0][kpoint] = identity;

        const std::vector<cplxdb> phases{
            std::polar(1.0, 0.2 * (kpoint + 1)),
            std::polar(1.0, -0.3 * (kpoint + 1))};
        if (kpoint % size_global == myid_global)
        {
            ComplexMatrix phased = identity;
            for (int band = 0; band < states; ++band)
            {
                for (int ao = 0; ao < states; ++ao)
                {
                    phased(band, ao) *=
                        phases[static_cast<std::size_t>(band)];
                }
            }
            distributed_basis.get_eigenvectors()[0][0][kpoint] =
                std::move(phased);
        }

        for (int direction = 0; direction < 3; ++direction)
        {
            auto& matrix = expected_velocity[0][kpoint][direction];
            matrix(0, 0) = 0.5 + kpoint + direction;
            matrix(1, 1) = -0.25 - 0.1 * kpoint;
            matrix(0, 1) = cplxdb(
                0.2 + 0.01 * kpoint, 0.3 + 0.02 * direction);
            matrix(1, 0) = std::conj(matrix(0, 1));
        }
    }

    source_velocity = expected_velocity;
    for (int kpoint = 0; kpoint < kpoints; ++kpoint)
    {
        const std::vector<cplxdb> phases{
            std::polar(1.0, 0.2 * (kpoint + 1)),
            std::polar(1.0, -0.3 * (kpoint + 1))};
        for (int direction = 0; direction < 3; ++direction)
        {
            for (int row = 0; row < states; ++row)
            {
                for (int column = 0; column < states; ++column)
                {
                    source_velocity[0][kpoint][direction](row, column) =
                        std::conj(phases[static_cast<std::size_t>(row)]) *
                        expected_velocity[0][kpoint][direction](row, column) *
                        phases[static_cast<std::size_t>(column)];
                }
            }
        }
    }

    const auto alignment = align_distributed_velocity_to_reference_wfc(
        distributed_basis, reference, source_velocity,
        mpi_comm_global_h);
    assert(alignment.maximum_relative_wfc_residual < 1.0e-14);
    assert(alignment.maximum_unitarity_residual < 1.0e-14);
    assert(alignment.maximum_transform_deviation_from_identity > 0.1);
    assert_velocity_equal(source_velocity, expected_velocity);

    MeanField missing_basis = distributed_basis;
    if (myid_global == 0)
    {
        missing_basis.get_eigenvectors()[0][0].erase(0);
    }
    VelocityMatrix untouched = expected_velocity;
    bool threw = false;
    try
    {
        (void)align_distributed_velocity_to_reference_wfc(
            missing_basis, reference, untouched,
            mpi_comm_global_h);
    }
    catch (const std::invalid_argument&)
    {
        threw = true;
    }
    assert(threw);
    assert_velocity_equal(untouched, expected_velocity);
}

} // namespace

int main(int argc, char* argv[])
{
    int provided = 0;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    init_global_mpi(MPI_COMM_WORLD);
    init_global_io();

    test_distributed_velocity_basis_alignment();

    finalize_global_io();
    finalize_global_mpi();
    MPI_Finalize();
    return 0;
}
