#include "../math/complexmatrix.h"
#include "testutils.h"
#include <random>
#include <ctime>

using namespace librpa_int;

void test_rectangular_complex_matrix_product()
{
    ComplexMatrix a(2, 3);
    ComplexMatrix b(3, 4);
    for (int i = 0; i < a.nr; ++i)
        for (int j = 0; j < a.nc; ++j)
            a(i, j) = {0.25 + 0.5 * i + 0.125 * j, -0.2 + 0.3 * i - 0.1 * j};
    for (int i = 0; i < b.nr; ++i)
        for (int j = 0; j < b.nc; ++j)
            b(i, j) = {-0.4 + 0.2 * i - 0.05 * j, 0.15 - 0.07 * i + 0.11 * j};

    const auto product = a * b;
    for (int i = 0; i < product.nr; ++i)
    {
        for (int j = 0; j < product.nc; ++j)
        {
            std::complex<double> expected{0.0, 0.0};
            for (int k = 0; k < a.nc; ++k) expected += a(i, k) * b(k, j);
            assert(std::abs(product(i, j) - expected) < 1e-12);
        }
    }
}

void test_power_hemat()
{
    ComplexMatrix a(2, 2);
    a(0, 0) = cmplx(1., 0.);
    a(1, 1) = cmplx(-1., 0.);
    // print_complex_matrix("a", a);
    auto b = power_hemat(a, 2, true);
    // print_complex_matrix("sq_a", b);
    assert( fequal(b(0, 0), cmplx(1., 0.)) &&
            fequal(b(1, 1), cmplx(1., 0.))
            );
    // check eigenvectors, with eigenvalue from big to small
    // print_complex_matrix("a", a);
    assert( fequal(a(0, 0), cmplx(1., 0.)) &&
            fequal(a(1, 1), cmplx(1., 0.))
            );

    a.zero_out();
    a(0, 0) = cmplx(-1., 0.);
    a(1, 1) = cmplx(1., 0.);
    power_hemat(a, 2, true);
    // check eigenvectors, with eigenvalue from big to small
    // print_complex_matrix("sq_a", b);
    assert(fequal(a(0, 1), cmplx(1., 0.)) &&
           fequal(a(1, 0), cmplx(1., 0.))
           );

    a.zero_out();
    a(0, 0) = cmplx(2., 0.);
    a(1, 1) = cmplx(2., 0.);
    a(0, 1) = cmplx(0., -2.);
    a(1, 0) = cmplx(0., 2.);
    // print_complex_matrix("a", a);
    b = power_hemat(a, 0.5, false);
    // print_complex_matrix("sqrt_a", b);
    // phase is undetermined, use square of sqrt to check
    /* assert( fequal(b(0, 0), cmplx(1, 0)) && */
    /*         fequal(b(1, 1), cmplx(1, 0)) && */
    /*         fequal(b(0, 1), cmplx(0, -1)) && */
    /*         fequal(b(1, 0), cmplx(0, 1)) */
    /*         ); */
    assert(is_matmul_AB_equal_C(2, 2, 2, b.c, b.c, a.c));

    ComplexMatrix iden(10, 10);
    iden.set_as_identity_matrix();
    // reset and randomize.
    a.create(10, 10);
    std::default_random_engine e(time(0));
    std::uniform_real_distribution<double> d(-1.0, 1.0);
    for (int i = 0; i < 10; i++)
    {
        a(i, i) = d(e);
        for (int j = i + 1; j < 10; j++)
        {
            a(i, j) = {d(e), d(e)};
            a(j, i) = conj(a(i, j));
        }
    }
    power_hemat(a, 2, true);
    // confirm a is unitary after solving and keep_ev set True
    const ComplexMatrix aconja = a * transpose(a, true);
    is_mat_A_equal_B(10, 10, aconja.c, iden.c, false, print_mat_equal, {1e-14, 0.0});
}

int main (int argc, char *argv[])
{
    test_rectangular_complex_matrix_product();
    test_power_hemat();
    
    return 0;
}
