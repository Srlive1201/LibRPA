#include "../complexmatrix.h"
#include "testutils.h"

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
    // print_complex_matrix("a", a);
    // check eigenvectors, with eigenvalue from small to big
    assert(fequal(a(0, 1), cmplx(1., 0.)) &&
           fequal(a(1, 0), cmplx(1., 0.))
           );

    a.zero_out();
    a(0, 0) = cmplx(-1., 0.);
    a(1, 1) = cmplx(1., 0.);
    power_hemat(a, 2, true);
    // print_complex_matrix("sq_a", b);
    assert( fequal(a(0, 0), cmplx(1., 0.)) &&
            fequal(a(1, 1), cmplx(1., 0.))
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
}

int main (int argc, char *argv[])
{
    test_power_hemat();
    
    return 0;
}
