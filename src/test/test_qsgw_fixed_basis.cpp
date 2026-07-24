#ifdef NDEBUG
#undef NDEBUG
#endif

#include "../qsgw/fixed_basis.h"

#include <cassert>
#include <cmath>
#include <complex>
#include <iostream>
#include <limits>
#include <stdexcept>

using librpa_int::ComplexMatrix;
using librpa_int::MeanField;
using librpa_int::Matz;
using librpa_int::cplxdb;
using librpa_int::qsgw::SpinKMatrixMap;
using librpa_int::qsgw::ScopedReferenceEigenvectors;
using librpa_int::qsgw::VelocityMatrix;
using librpa_int::qsgw::align_velocity_to_reference_wfc;
using librpa_int::qsgw::diagonalize_in_reference_basis;
using librpa_int::qsgw::prepare_fhi_aims_interband_velocity;

namespace
{

void assert_close(const cplxdb actual, const cplxdb expected,
                  const double tolerance = 1.0e-12)
{
    assert(std::abs(actual - expected) < tolerance);
}

void assert_meanfield_equal(const MeanField& actual, const MeanField& expected)
{
    assert(actual.get_n_spins() == expected.get_n_spins());
    assert(actual.get_n_kpoints() == expected.get_n_kpoints());
    assert(actual.get_n_bands() == expected.get_n_bands());
    assert(actual.get_n_aos() == expected.get_n_aos());
    assert(actual.get_n_spinor() == expected.get_n_spinor());
    assert_close(actual.get_efermi(), expected.get_efermi());

    for (int spin = 0; spin < actual.get_n_spins(); ++spin)
    {
        for (int kpoint = 0; kpoint < actual.get_n_kpoints(); ++kpoint)
        {
            for (int band = 0; band < actual.get_n_bands(); ++band)
            {
                assert_close(actual.get_eigenvals()[spin](kpoint, band),
                             expected.get_eigenvals()[spin](kpoint, band));
                assert_close(actual.get_weight()[spin](kpoint, band),
                             expected.get_weight()[spin](kpoint, band));
            }
            for (int spinor = 0; spinor < actual.get_n_spinor(); ++spinor)
            {
                const auto& actual_wfc =
                    actual.get_eigenvectors().at(spin).at(spinor).at(kpoint);
                const auto& expected_wfc =
                    expected.get_eigenvectors().at(spin).at(spinor).at(kpoint);
                assert(actual_wfc.nr == expected_wfc.nr);
                assert(actual_wfc.nc == expected_wfc.nc);
                for (int row = 0; row < actual_wfc.nr; ++row)
                {
                    for (int column = 0; column < actual_wfc.nc; ++column)
                    {
                        assert_close(actual_wfc(row, column),
                                     expected_wfc(row, column));
                    }
                }
            }
        }
    }
}

void assert_velocity_equal(const VelocityMatrix& actual,
                           const VelocityMatrix& expected)
{
    assert(actual.size() == expected.size());
    for (std::size_t spin = 0; spin < actual.size(); ++spin)
    {
        assert(actual[spin].size() == expected[spin].size());
        for (std::size_t kpoint = 0; kpoint < actual[spin].size(); ++kpoint)
        {
            assert(actual[spin][kpoint].size() == expected[spin][kpoint].size());
            for (std::size_t direction = 0;
                 direction < actual[spin][kpoint].size(); ++direction)
            {
                const auto& actual_matrix = actual[spin][kpoint][direction];
                const auto& expected_matrix = expected[spin][kpoint][direction];
                assert(actual_matrix.nr == expected_matrix.nr);
                assert(actual_matrix.nc == expected_matrix.nc);
                for (int row = 0; row < actual_matrix.nr; ++row)
                {
                    for (int column = 0; column < actual_matrix.nc; ++column)
                    {
                        assert_close(actual_matrix(row, column),
                                     expected_matrix(row, column));
                    }
                }
            }
        }
    }
}

void initialize_velocity_matrix(VelocityMatrix& velocity,
                                const int spins,
                                const int kpoints,
                                const int states)
{
    velocity.assign(spins, {});
    for (int spin = 0; spin < spins; ++spin)
    {
        velocity[spin].assign(kpoints, {});
        for (int kpoint = 0; kpoint < kpoints; ++kpoint)
        {
            velocity[spin][kpoint].assign(3, ComplexMatrix(states, states));
        }
    }
}

template <typename Function>
void assert_throws(Function&& function)
{
    bool threw = false;
    try
    {
        function();
    }
    catch (const std::exception&)
    {
        threw = true;
    }
    assert(threw);
}

Matz collect_wfc_rows(const MeanField& meanfield)
{
    Matz output(meanfield.get_n_bands(),
                meanfield.get_n_aos() * meanfield.get_n_spinor());
    for (int band = 0; band < meanfield.get_n_bands(); ++band)
    {
        for (int spinor = 0; spinor < meanfield.get_n_spinor(); ++spinor)
        {
            const auto& block = meanfield.get_eigenvectors().at(0).at(spinor).at(0);
            for (int ao = 0; ao < meanfield.get_n_aos(); ++ao)
            {
                output(band, ao * meanfield.get_n_spinor() + spinor) =
                    block(band, ao);
            }
        }
    }
    return output;
}

void initialize_identity_wfc(MeanField& meanfield, const int kpoint = 0)
{
    ComplexMatrix identity(2, 2);
    identity(0, 0) = 1.0;
    identity(1, 1) = 1.0;
    meanfield.get_eigenvectors()[0][0][kpoint] = identity;
}

void test_scoped_reference_eigenvectors_restores_live_state()
{
    MeanField reference(1, 1, 2, 2, 1);
    initialize_identity_wfc(reference);
    reference.get_eigenvals()[0](0, 0) = -0.5;
    reference.get_eigenvals()[0](0, 1) = 0.5;
    const MeanField reference_before = reference;

    MeanField live = reference;
    live.get_eigenvals()[0](0, 0) = -0.25;
    live.get_eigenvals()[0](0, 1) = 0.75;
    live.get_weight()[0](0, 0) = 2.0;
    live.get_efermi() = 0.125;
    auto& live_wfc = live.get_eigenvectors()[0][0][0];
    live_wfc(0, 0) = 0.0;
    live_wfc(0, 1) = 1.0;
    live_wfc(1, 0) = -1.0;
    live_wfc(1, 1) = 0.0;
    const MeanField live_before = live;

    bool threw = false;
    try
    {
        ScopedReferenceEigenvectors scope(live, reference);
        const Matz projected_wfc = collect_wfc_rows(live);
        const Matz reference_wfc = collect_wfc_rows(reference);
        for (int row = 0; row < 2; ++row)
        {
            for (int column = 0; column < 2; ++column)
            {
                assert_close(projected_wfc(row, column),
                             reference_wfc(row, column));
            }
        }
        assert_close(live.get_eigenvals()[0](0, 0), -0.25);
        assert_close(live.get_eigenvals()[0](0, 1), 0.75);
        assert_close(live.get_weight()[0](0, 0), 2.0);
        assert_close(live.get_efermi(), 0.125);
        throw std::runtime_error("exercise exception restoration");
    }
    catch (const std::runtime_error&)
    {
        threw = true;
    }
    assert(threw);
    assert_meanfield_equal(live, live_before);
    assert_meanfield_equal(reference, reference_before);

    assert_throws([&] {
        ScopedReferenceEigenvectors invalid(reference, reference);
    });

    MeanField mismatched(1, 1, 1, 1, 1);
    mismatched.get_eigenvectors()[0][0][0] = ComplexMatrix(1, 1);
    const MeanField live_before_mismatch = live;
    assert_throws([&] {
        ScopedReferenceEigenvectors invalid(live, mismatched);
    });
    assert_meanfield_equal(live, live_before_mismatch);
}

void test_fixed_reference_updates_live_wfc_from_mf0_each_iteration()
{
    MeanField reference(1, 1, 2, 2, 1);
    initialize_identity_wfc(reference);
    reference.get_eigenvals()[0](0, 0) = -0.5;
    reference.get_eigenvals()[0](0, 1) = 0.5;
    reference.get_weight()[0](0, 0) = 2.0;
    reference.get_efermi() = 0.25;
    const MeanField reference_before = reference;

    MeanField live = reference;
    SpinKMatrixMap hamiltonian;
    hamiltonian[0][0] = Matz(2, 2);
    hamiltonian[0][0](0, 0) = 1.0;
    hamiltonian[0][0](0, 1) = cplxdb(0.0, 0.4);
    hamiltonian[0][0](1, 0) = cplxdb(0.0, -0.4);
    hamiltonian[0][0](1, 1) = 2.0;

    const auto result = diagonalize_in_reference_basis(
        live, reference, hamiltonian);
    const Matz& unitary = result.unitary.at(0).at(0);
    const Matz unitary_transpose = transpose(unitary, false);
    const Matz unitary_dagger = transpose(unitary, true);

    // MeanField stores ket coefficients as band-by-AO rows, so a basis
    // rotation whose eigenvectors are columns of U is C_live = U^T C_0.
    const Matz expected_wfc =
        unitary_transpose * collect_wfc_rows(reference);
    const Matz actual_wfc = collect_wfc_rows(live);
    for (int row = 0; row < 2; ++row)
    {
        for (int column = 0; column < 2; ++column)
        {
            assert_close(actual_wfc(row, column), expected_wfc(row, column));
        }
    }

    assert_meanfield_equal(reference, reference_before);

    // The immutable reference remains the producer basis after live updates.
    const Matz reference_after_wfc = collect_wfc_rows(reference);
    const Matz reference_before_wfc = collect_wfc_rows(reference_before);
    for (int row = 0; row < 2; ++row)
    {
        for (int column = 0; column < 2; ++column)
        {
            assert_close(reference_after_wfc(row, column),
                         reference_before_wfc(row, column));
        }
    }

    const Matz diagonal = unitary_dagger * hamiltonian.at(0).at(0) * unitary;
    assert_close(diagonal(0, 1), 0.0);
    assert_close(diagonal(1, 0), 0.0);
    assert_close(live.get_eigenvals()[0](0, 0), diagonal(0, 0).real());
    assert_close(live.get_eigenvals()[0](0, 1), diagonal(1, 1).real());

    // A second iteration must rotate directly from the immutable reference,
    // not compound the first iteration's live wavefunction.
    hamiltonian[0][0](0, 0) = 2.0;
    hamiltonian[0][0](0, 1) = cplxdb(0.3, 0.2);
    hamiltonian[0][0](1, 0) = cplxdb(0.3, -0.2);
    hamiltonian[0][0](1, 1) = 0.5;

    const auto second_result = diagonalize_in_reference_basis(
        live, reference, hamiltonian);
    const Matz& second_unitary = second_result.unitary.at(0).at(0);
    const Matz second_unitary_transpose =
        transpose(second_unitary, false);
    const Matz second_unitary_dagger = transpose(second_unitary, true);
    const Matz second_expected_wfc =
        second_unitary_transpose * collect_wfc_rows(reference);
    const Matz second_actual_wfc = collect_wfc_rows(live);
    for (int row = 0; row < 2; ++row)
    {
        for (int column = 0; column < 2; ++column)
        {
            assert_close(second_actual_wfc(row, column),
                         second_expected_wfc(row, column));
        }
    }

    const Matz second_diagonal =
        second_unitary_dagger * hamiltonian.at(0).at(0) * second_unitary;
    assert_close(second_diagonal(0, 1), 0.0);
    assert_close(second_diagonal(1, 0), 0.0);
    assert_close(live.get_eigenvals()[0](0, 0),
                 second_diagonal(0, 0).real());
    assert_close(live.get_eigenvals()[0](0, 1),
                 second_diagonal(1, 1).real());

    const Matz reference_after_second = collect_wfc_rows(reference);
    for (int row = 0; row < 2; ++row)
    {
        for (int column = 0; column < 2; ++column)
        {
            assert_close(reference_after_second(row, column),
                         reference_before_wfc(row, column));
        }
    }
    assert_meanfield_equal(reference, reference_before);
    assert_close(live.get_weight()[0](0, 0), 2.0);
    assert_close(live.get_weight()[0](0, 1), 0.0);
    assert_close(live.get_efermi(), 0.25);
}

void test_invalid_late_kpoint_does_not_partially_update_live_state()
{
    MeanField reference(1, 2, 2, 2, 1);
    initialize_identity_wfc(reference, 0);
    initialize_identity_wfc(reference, 1);
    MeanField live = reference;
    live.get_eigenvals()[0](0, 0) = -1.0;
    live.get_eigenvals()[0](0, 1) = 1.0;
    live.get_eigenvals()[0](1, 0) = -2.0;
    live.get_eigenvals()[0](1, 1) = 2.0;
    const MeanField live_before = live;

    SpinKMatrixMap hamiltonian;
    hamiltonian[0][0] = Matz(2, 2);
    hamiltonian[0][0](0, 0) = 0.5;
    hamiltonian[0][0](1, 1) = 1.5;
    hamiltonian[0][1] = Matz(1, 1);
    hamiltonian[0][1](0, 0) = 3.0;

    assert_throws([&] {
        (void)diagonalize_in_reference_basis(
            live, reference, hamiltonian);
    });

    assert_meanfield_equal(live, live_before);

    SpinKMatrixMap valid_hamiltonian = hamiltonian;
    valid_hamiltonian[0][1] = Matz(2, 2);
    valid_hamiltonian[0][1](0, 0) = 2.0;
    valid_hamiltonian[0][1](1, 1) = 3.0;
    assert_throws([&] {
        (void)diagonalize_in_reference_basis(
            reference, reference, valid_hamiltonian);
    });
}

void test_invalid_matrix_data_does_not_mutate_live_state()
{
    MeanField reference(1, 1, 2, 2, 1);
    initialize_identity_wfc(reference);
    reference.get_eigenvals()[0](0, 0) = -0.5;
    reference.get_eigenvals()[0](0, 1) = 0.5;
    MeanField live = reference;

    const auto assert_rejected_without_mutation = [&](const Matz& matrix) {
        SpinKMatrixMap hamiltonian;
        hamiltonian[0][0] = matrix;
        const MeanField live_before = live;
        assert_throws([&] {
            (void)diagonalize_in_reference_basis(
                live, reference, hamiltonian);
        });
        assert_meanfield_equal(live, live_before);
    };

    Matz nonhermitian(2, 2);
    nonhermitian(0, 0) = 1.0;
    nonhermitian(0, 1) = cplxdb(0.0, 0.4);
    nonhermitian(1, 0) = cplxdb(0.0, 0.4);
    nonhermitian(1, 1) = 2.0;
    assert_rejected_without_mutation(nonhermitian);

    Matz nonfinite(2, 2);
    nonfinite(0, 0) = std::numeric_limits<double>::quiet_NaN();
    nonfinite(1, 1) = 2.0;
    assert_rejected_without_mutation(nonfinite);
}

void test_velocity_basis_unitary_is_aligned_to_fixed_reference()
{
    MeanField reference(1, 1, 2, 2, 1);
    initialize_identity_wfc(reference);
    MeanField source_basis = reference;
    const std::vector<cplxdb> phases{
        cplxdb(0.0, 1.0), cplxdb(0.6, -0.8)};
    for (int band = 0; band < 2; ++band)
    {
        for (int ao = 0; ao < 2; ++ao)
        {
            source_basis.get_eigenvectors()[0][0][0](band, ao) *=
                phases[static_cast<std::size_t>(band)];
        }
    }

    VelocityMatrix reference_velocity;
    initialize_velocity_matrix(reference_velocity, 1, 1, 2);
    reference_velocity[0][0][0](0, 0) = 1.0;
    reference_velocity[0][0][0](0, 1) = cplxdb(0.2, 0.3);
    reference_velocity[0][0][0](1, 0) = cplxdb(0.2, -0.3);
    reference_velocity[0][0][0](1, 1) = -0.5;
    reference_velocity[0][0][1](0, 1) = cplxdb(-0.4, 0.1);
    reference_velocity[0][0][1](1, 0) = cplxdb(-0.4, -0.1);
    reference_velocity[0][0][2](0, 0) = 0.75;
    reference_velocity[0][0][2](1, 1) = 1.25;

    VelocityMatrix source_velocity = reference_velocity;
    for (int direction = 0; direction < 3; ++direction)
    {
        for (int row = 0; row < 2; ++row)
        {
            for (int column = 0; column < 2; ++column)
            {
                source_velocity[0][0][direction](row, column) =
                    std::conj(phases[static_cast<std::size_t>(row)]) *
                    reference_velocity[0][0][direction](row, column) *
                    phases[static_cast<std::size_t>(column)];
            }
        }
    }

    const auto alignment = align_velocity_to_reference_wfc(
        source_basis, reference, source_velocity);
    assert(alignment.maximum_relative_wfc_residual < 1.0e-14);
    assert(alignment.maximum_unitarity_residual < 1.0e-14);
    assert(alignment.maximum_basis_inverse_residual < 1.0e-14);
    assert(alignment.maximum_transform_deviation_from_identity > 1.0);
    assert_velocity_equal(source_velocity, reference_velocity);

    MeanField mixed_basis = reference;
    const double inverse_sqrt_two = 1.0 / std::sqrt(2.0);
    mixed_basis.get_eigenvectors()[0][0][0](0, 0) = inverse_sqrt_two;
    mixed_basis.get_eigenvectors()[0][0][0](0, 1) = inverse_sqrt_two;
    mixed_basis.get_eigenvectors()[0][0][0](1, 0) = -inverse_sqrt_two;
    mixed_basis.get_eigenvectors()[0][0][0](1, 1) = inverse_sqrt_two;

    ComplexMatrix transform(2, 2);
    transform(0, 0) = inverse_sqrt_two;
    transform(0, 1) = inverse_sqrt_two;
    transform(1, 0) = -inverse_sqrt_two;
    transform(1, 1) = inverse_sqrt_two;
    VelocityMatrix mixed_velocity = reference_velocity;
    for (int direction = 0; direction < 3; ++direction)
    {
        mixed_velocity[0][0][direction] =
            librpa_int::conj(transform) *
            reference_velocity[0][0][direction] *
            transpose(transform, false);
    }
    const auto mixed_alignment = align_velocity_to_reference_wfc(
        mixed_basis, reference, mixed_velocity);
    assert(mixed_alignment.maximum_relative_wfc_residual < 1.0e-14);
    assert(mixed_alignment.maximum_unitarity_residual < 1.0e-14);
    assert(mixed_alignment.maximum_transform_deviation_from_identity > 0.7);
    assert_velocity_equal(mixed_velocity, reference_velocity);

    MeanField rounded_basis = reference;
    rounded_basis.get_eigenvectors()[0][0][0](0, 0) = 1.0 + 2.0e-10;
    rounded_basis.get_eigenvectors()[0][0][0](1, 1) = 1.0 - 1.0e-10;
    VelocityMatrix rounded_velocity = reference_velocity;
    const auto rounded_alignment = align_velocity_to_reference_wfc(
        rounded_basis, reference, rounded_velocity);
    assert(rounded_alignment.maximum_relative_wfc_residual > 1.0e-11);
    assert(rounded_alignment.maximum_relative_wfc_residual < 1.0e-8);
    assert(rounded_alignment.maximum_unitarity_residual < 1.0e-14);
    assert(rounded_alignment.maximum_raw_unitarity_residual > 1.0e-10);
    assert(rounded_alignment.maximum_unitary_projection_correction > 1.0e-11);
    assert_velocity_equal(rounded_velocity, reference_velocity);

    MeanField nonunitary_basis = reference;
    nonunitary_basis.get_eigenvectors()[0][0][0](0, 0) = 2.0;
    VelocityMatrix rejected_velocity = reference_velocity;
    const VelocityMatrix before_rejection = rejected_velocity;
    assert_throws([&] {
        (void)align_velocity_to_reference_wfc(
            nonunitary_basis, reference, rejected_velocity);
    });
    assert_velocity_equal(rejected_velocity, before_rejection);
}

void test_fhi_aims_interband_velocity_is_prepared_in_qsgw_only()
{
    MeanField reference(1, 1, 2, 2, 1);
    VelocityMatrix velocity;
    initialize_velocity_matrix(velocity, 1, 1, 2);
    for (int direction = 0; direction < 3; ++direction)
    {
        const cplxdb off_diagonal(
            0.25 + direction, -0.5 + 0.1 * direction);
        velocity[0][0][direction](0, 0) =
            cplxdb(1.0 + direction, 0.4);
        velocity[0][0][direction](0, 1) = off_diagonal;
        velocity[0][0][direction](1, 0) = std::conj(off_diagonal);
        velocity[0][0][direction](1, 1) =
            cplxdb(-2.0 - direction, -0.7);
    }
    const VelocityMatrix before = velocity;

    prepare_fhi_aims_interband_velocity(velocity, reference);
    for (int direction = 0; direction < 3; ++direction)
    {
        assert_close(velocity[0][0][direction](0, 0), 0.0);
        assert_close(velocity[0][0][direction](1, 1), 0.0);
        assert_close(velocity[0][0][direction](0, 1),
                     before[0][0][direction](0, 1));
        assert_close(velocity[0][0][direction](1, 0),
                     before[0][0][direction](1, 0));
    }

    VelocityMatrix invalid = before;
    invalid[0][0][2](1, 0) += cplxdb(0.0, 0.25);
    const VelocityMatrix invalid_before = invalid;
    assert_throws([&] {
        prepare_fhi_aims_interband_velocity(invalid, reference);
    });
    assert_velocity_equal(invalid, invalid_before);
}

} // namespace

int main()
{
    test_scoped_reference_eigenvectors_restores_live_state();
    test_fixed_reference_updates_live_wfc_from_mf0_each_iteration();
    test_invalid_late_kpoint_does_not_partially_update_live_state();
    test_invalid_matrix_data_does_not_mutate_live_state();
    test_velocity_basis_unitary_is_aligned_to_fixed_reference();
    test_fhi_aims_interband_velocity_is_prepared_in_qsgw_only();
    std::cout << "test_qsgw_fixed_basis: all tests passed\n";
    return 0;
}
