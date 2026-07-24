#ifdef NDEBUG
#undef NDEBUG
#endif

#include "../qsgw/effective_hamiltonian.h"

#include <cassert>
#include <complex>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

using librpa_int::Matz;
using librpa_int::cplxdb;
using librpa_int::qsgw::SpinKMatrixMap;
using librpa_int::qsgw::assemble_effective_hamiltonian;
using librpa_int::qsgw::build_reference_hamiltonian;

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

template <typename Function>
std::string exception_message(Function&& function)
{
    try
    {
        function();
    }
    catch (const std::exception& error)
    {
        return error.what();
    }
    assert(false);
    return {};
}

Matz hermitian(const double diagonal0, const double diagonal1,
               const cplxdb off_diagonal,
               const librpa_int::MAJOR major = librpa_int::MAJOR::ROW)
{
    Matz result(2, 2, major);
    result(0, 0) = diagonal0;
    result(1, 1) = diagonal1;
    result(0, 1) = off_diagonal;
    result(1, 0) = std::conj(off_diagonal);
    return result;
}

SpinKMatrixMap one_block(const Matz& value)
{
    SpinKMatrixMap result;
    result[0][0] = value.copy();
    return result;
}

void assert_close(const cplxdb actual, const cplxdb expected,
                  const double tolerance = 1.0e-13)
{
    assert(std::abs(actual - expected) < tolerance);
}

void test_qsgw_formula()
{
    const auto hks = one_block(hermitian(
        1.0, 2.0, {0.1, -0.2}, librpa_int::MAJOR::ROW));
    const auto vxc = one_block(hermitian(
        0.3, 0.4, {-0.2, 0.1}, librpa_int::MAJOR::COL));
    const auto exchange = one_block(hermitian(
        -0.5, -0.6, {0.05, 0.3}, librpa_int::MAJOR::COL));
    const auto correlation = one_block(hermitian(
        0.2, 0.1, {0.15, -0.1}, librpa_int::MAJOR::ROW));

    const Matz hks_before = hks.at(0).at(0).copy();
    const Matz vxc_before = vxc.at(0).at(0).copy();
    const Matz exchange_before = exchange.at(0).at(0).copy();
    const Matz correlation_before = correlation.at(0).at(0).copy();
    Matz expected(2, 2);
    for (int row = 0; row < 2; ++row)
    {
        for (int column = 0; column < 2; ++column)
        {
            expected(row, column) =
                hks.at(0).at(0)(row, column) -
                vxc.at(0).at(0)(row, column) +
                exchange.at(0).at(0)(row, column) +
                correlation.at(0).at(0)(row, column);
        }
    }

    const auto result = assemble_effective_hamiltonian(
        hks, vxc, exchange, correlation);
    const auto& matrix = result.at(0).at(0);
    for (int row = 0; row < 2; ++row)
    {
        for (int column = 0; column < 2; ++column)
        {
            assert_close(matrix(row, column), expected(row, column));
            assert_close(hks.at(0).at(0)(row, column),
                         hks_before(row, column));
            assert_close(vxc.at(0).at(0)(row, column),
                         vxc_before(row, column));
            assert_close(exchange.at(0).at(0)(row, column),
                         exchange_before(row, column));
            assert_close(correlation.at(0).at(0)(row, column),
                         correlation_before(row, column));
        }
    }
}

void test_layout_and_finiteness_mismatches_are_rejected()
{
    const auto valid = one_block(hermitian(1.0, 2.0, {0.0, 0.0}));
    auto missing = valid;
    missing.at(0).erase(0);
    assert_throws([&] {
        assemble_effective_hamiltonian(
            valid, missing, valid, valid);
    });

    auto nonfinite = one_block(valid.at(0).at(0));
    nonfinite.at(0).at(0)(0, 0) =
        std::numeric_limits<double>::quiet_NaN();
    assert(std::isfinite(valid.at(0).at(0)(0, 0).real()));
    assert_throws([&] {
        assemble_effective_hamiltonian(
            valid, valid, valid, nonfinite);
    });
}

void test_nonhermitian_components_use_the_legacy_upper_triangle()
{
    const auto valid = one_block(hermitian(1.0, 2.0, {0.0, 0.0}));
    auto exchange = one_block(hermitian(-0.5, -0.6, {0.05, 0.3}));
    exchange.at(0).at(0)(0, 0) += cplxdb(0.0, 4.0e-7);
    exchange.at(0).at(0)(1, 0) += cplxdb(0.0, 2.5e-7);

    const auto result = assemble_effective_hamiltonian(
        valid, valid, exchange, valid);
    const auto& matrix = result.at(0).at(0);
    assert_close(matrix(0, 0), cplxdb(0.5, 0.0));
    assert_close(matrix(0, 1), cplxdb(0.05, 0.3));
    assert_close(matrix(1, 0), std::conj(matrix(0, 1)));
}

void test_reference_hamiltonian_is_the_initial_eigenvalue_diagonal()
{
    librpa_int::MeanField reference(1, 2, 3, 1, 1);
    for (int kpoint = 0; kpoint < 2; ++kpoint)
    {
        for (int band = 0; band < 3; ++band)
        {
            reference.get_eigenvals()[0](kpoint, band) =
                10.0 * kpoint + band + 0.25;
        }
    }
    const auto result = build_reference_hamiltonian(reference);
    for (int kpoint = 0; kpoint < 2; ++kpoint)
    {
        for (int row = 0; row < 3; ++row)
        {
            for (int column = 0; column < 3; ++column)
            {
                const cplxdb expected = row == column
                    ? cplxdb(10.0 * kpoint + row + 0.25, 0.0)
                    : cplxdb(0.0, 0.0);
                assert_close(result.at(0).at(kpoint)(row, column), expected);
            }
        }
    }
}

} // namespace

int main()
{
    test_qsgw_formula();
    test_layout_and_finiteness_mismatches_are_rejected();
    test_nonhermitian_components_use_the_legacy_upper_triangle();
    test_reference_hamiltonian_is_the_initial_eigenvalue_diagonal();
    std::cout << "test_qsgw_effective_hamiltonian: all tests passed\n";
    return 0;
}
