#ifdef NDEBUG
#undef NDEBUG
#endif

#include "../qsgw/hamiltonian_cut.h"

#include <cassert>
#include <cmath>
#include <complex>
#include <iostream>
#include <limits>
#include <stdexcept>

using librpa_int::MAJOR;
using librpa_int::Matz;
using librpa_int::MeanField;
using librpa_int::cplxdb;
using librpa_int::qsgw::HamiltonianCutMode;
using librpa_int::qsgw::HamiltonianCutOptions;
using librpa_int::qsgw::SpinKMatrixMap;
using librpa_int::qsgw::apply_hamiltonian_cut;
using librpa_int::qsgw::hamiltonian_cut_mode_from_int;

namespace
{

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

void assert_close(const cplxdb actual, const cplxdb expected,
                  const double tolerance = 1.0e-13)
{
    assert(std::abs(actual - expected) < tolerance);
}

Matz make_matrix(const int dimension, const double diagonal_offset)
{
    Matz result(dimension, dimension, MAJOR::ROW);
    for (int row = 0; row < dimension; ++row)
    {
        result(row, row) = diagonal_offset + row;
        for (int column = row + 1; column < dimension; ++column)
        {
            const cplxdb value(0.1 * (row + 1), 0.05 * (column + 1));
            result(row, column) = value;
            result(column, row) = std::conj(value);
        }
    }
    return result;
}

MeanField make_live_meanfield()
{
    MeanField live(1, 1, 5, 1, 1);
    live.get_efermi() = 0.0;
    live.get_eigenvals()[0](0, 0) = -2.0;
    live.get_eigenvals()[0](0, 1) = -0.5;
    live.get_eigenvals()[0](0, 2) = 0.0;
    live.get_eigenvals()[0](0, 3) = 0.4;
    live.get_eigenvals()[0](0, 4) = 0.8;
    return live;
}

void test_mode_zero_is_an_exact_uncut_copy()
{
    MeanField live = make_live_meanfield();
    live.get_efermi() = std::numeric_limits<double>::quiet_NaN();
    live.get_eigenvals()[0](0, 0) =
        std::numeric_limits<double>::quiet_NaN();
    SpinKMatrixMap raw;
    SpinKMatrixMap reference;
    raw[0][0] = make_matrix(5, 10.0);
    reference[0][0] = make_matrix(5, -3.0);
    const Matz raw_before = raw.at(0).at(0).copy();

    HamiltonianCutOptions options;
    options.mode = HamiltonianCutMode::Uncut;
    options.unoccupied_keep = 0;
    options.shift_ha = 123.0;
    const auto result = apply_hamiltonian_cut(raw, reference, live, options);

    for (int row = 0; row < 5; ++row)
    {
        for (int column = 0; column < 5; ++column)
        {
            assert_close(result.at(0).at(0)(row, column),
                         raw_before(row, column));
            assert_close(raw.at(0).at(0)(row, column),
                         raw_before(row, column));
        }
    }
}

void test_cut_modes_require_a_finite_live_spectrum()
{
    MeanField live = make_live_meanfield();
    SpinKMatrixMap raw;
    SpinKMatrixMap reference;
    raw[0][0] = make_matrix(5, 10.0);
    reference[0][0] = make_matrix(5, -3.0);

    HamiltonianCutOptions options;
    options.mode = HamiltonianCutMode::ReferenceDiagonal;
    live.get_efermi() = std::numeric_limits<double>::quiet_NaN();
    assert_throws([&] {
        (void)apply_hamiltonian_cut(raw, reference, live, options);
    });

    live = make_live_meanfield();
    live.get_eigenvals()[0](0, 2) =
        std::numeric_limits<double>::quiet_NaN();
    assert_throws([&] {
        (void)apply_hamiltonian_cut(raw, reference, live, options);
    });
}

void test_mode_one_restores_reference_above_legacy_active_limit()
{
    const MeanField live = make_live_meanfield();
    SpinKMatrixMap raw;
    SpinKMatrixMap reference;
    raw[0][0] = make_matrix(5, 10.0);
    reference[0][0] = make_matrix(5, -3.0);

    HamiltonianCutOptions options;
    options.mode = HamiltonianCutMode::ReferenceDiagonal;
    options.unoccupied_keep = 1;
    const auto result = apply_hamiltonian_cut(raw, reference, live, options);
    const Matz& cut = result.at(0).at(0);

    // Strict energy < efermi gives Nocc=2. Therefore bands 0,1,2 survive
    // and the cut starts at index 3, matching the frozen 8476213 source.
    for (int row = 0; row < 3; ++row)
        for (int column = 0; column < 3; ++column)
            assert_close(cut(row, column), raw.at(0).at(0)(row, column));
    for (int row = 0; row < 5; ++row)
    {
        for (int column = 0; column < 5; ++column)
        {
            if (row < 3 && column < 3) continue;
            const cplxdb expected = row == column
                ? reference.at(0).at(0)(row, row)
                : cplxdb(0.0, 0.0);
            assert_close(cut(row, column), expected);
        }
    }
}

void test_mode_two_adds_hartree_shift_to_reference_diagonal()
{
    const MeanField live = make_live_meanfield();
    SpinKMatrixMap raw;
    SpinKMatrixMap reference;
    raw[0][0] = make_matrix(5, 10.0);
    reference[0][0] = make_matrix(5, -3.0);

    HamiltonianCutOptions options;
    options.mode = HamiltonianCutMode::ShiftedReferenceDiagonal;
    options.unoccupied_keep = 1;
    options.shift_ha = 20.0;
    const auto result = apply_hamiltonian_cut(raw, reference, live, options);
    const Matz& cut = result.at(0).at(0);

    assert_close(cut(2, 2), raw.at(0).at(0)(2, 2));
    assert_close(cut(3, 3), reference.at(0).at(0)(3, 3) + 20.0);
    assert_close(cut(4, 4), reference.at(0).at(0)(4, 4) + 20.0);
    assert_close(cut(2, 3), 0.0);
    assert_close(cut(3, 2), 0.0);
}

void test_keep_covering_all_bands_is_uncut()
{
    const MeanField live = make_live_meanfield();
    SpinKMatrixMap raw;
    SpinKMatrixMap reference;
    raw[0][0] = make_matrix(5, 10.0);
    reference[0][0] = make_matrix(5, -3.0);

    HamiltonianCutOptions options;
    options.mode = HamiltonianCutMode::ShiftedReferenceDiagonal;
    options.unoccupied_keep = 5;
    const auto result = apply_hamiltonian_cut(raw, reference, live, options);
    for (int row = 0; row < 5; ++row)
        for (int column = 0; column < 5; ++column)
            assert_close(result.at(0).at(0)(row, column),
                         raw.at(0).at(0)(row, column));
}

void test_invalid_options_and_layouts_fail_closed()
{
    const MeanField live = make_live_meanfield();
    SpinKMatrixMap raw;
    SpinKMatrixMap reference;
    raw[0][0] = make_matrix(5, 10.0);
    reference[0][0] = make_matrix(5, -3.0);

    assert(hamiltonian_cut_mode_from_int(0) == HamiltonianCutMode::Uncut);
    assert(hamiltonian_cut_mode_from_int(1) ==
           HamiltonianCutMode::ReferenceDiagonal);
    assert(hamiltonian_cut_mode_from_int(2) ==
           HamiltonianCutMode::ShiftedReferenceDiagonal);
    assert_throws([&] { (void)hamiltonian_cut_mode_from_int(-1); });
    assert_throws([&] { (void)hamiltonian_cut_mode_from_int(3); });

    HamiltonianCutOptions options;
    options.unoccupied_keep = -1;
    assert_throws([&] {
        (void)apply_hamiltonian_cut(raw, reference, live, options);
    });
    options.unoccupied_keep = 1;
    options.shift_ha = std::numeric_limits<double>::infinity();
    assert_throws([&] {
        (void)apply_hamiltonian_cut(raw, reference, live, options);
    });
    options.shift_ha = 20.0;
    SpinKMatrixMap wrong_reference = reference;
    wrong_reference.at(0).at(0) = make_matrix(4, -3.0);
    assert_throws([&] {
        (void)apply_hamiltonian_cut(raw, wrong_reference, live, options);
    });
}

} // namespace

int main()
{
    test_mode_zero_is_an_exact_uncut_copy();
    test_cut_modes_require_a_finite_live_spectrum();
    test_mode_one_restores_reference_above_legacy_active_limit();
    test_mode_two_adds_hartree_shift_to_reference_diagonal();
    test_keep_covering_all_bands_is_uncut();
    test_invalid_options_and_layouts_fail_closed();
    std::cout << "test_qsgw_hamiltonian_cut: all tests passed\n";
    return 0;
}
