#include "../matrix.h"
#include <stdexcept>
#include <cstdlib>
#include <ctime>
#include <cmath>

bool double_equal(double a, double b)
{
    double eqthres = 1e-14;
    return fabs(a - b) < eqthres;
}

void test_power_symat()
{
    matrix a(2, 2);
    a(0, 0) = 2;
    a(1, 1) = 1;
    print_matrix("a", a);
    matrix b = power_symat(a, 2);
    print_matrix("sq_a", b);
    assert( double_equal(b(0, 0), 4) &&
            double_equal(b(1, 1), 1) );
    a = power_symat(b, 0.5);
    assert( double_equal(a(0, 0), 2) &&
            double_equal(a(1, 1), 1) );

    a.zero_out();
    a(0, 1) = 1;
    a(1, 0) = 1;
    print_matrix("a", a);
    b = power_symat(a, 2);
    print_matrix("sq_a", b);
    assert( double_equal(b(0, 0), 1) &&
            double_equal(b(1, 1), 1) &&
            double_equal(b(1, 0), 0) &&
            double_equal(b(0, 1), 0));
}

int main (int argc, char *argv[])
{
    using namespace std;
    matrix m1(5, 6, true);
    srand(time(0));
    for (int i = 0; i < m1.size; i++)
        m1.c[i] = rand();

    matrix m2(m1);
    assert(m1.size == m2.size);
    for (int i = 0; i < m1.size; i++)
        if ( m1.c[i] != m2.c[i] ) throw logic_error("copy constructor fail");
    if (! (m1 == m2 )) throw logic_error("direct compare fail");

    test_power_symat();

    return 0;

}

