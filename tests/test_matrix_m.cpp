#include "../src/matrix_m.h"
#include "testutils.h"

void test_constuctor()
{
    matrix_m<double> m1(2, 3, MAJOR::COL);
    m1.zero_out();
    m1(0, 1) = 1.0;
    m1(1, 0) = 2.0;
    printf("Matrix (set column-major)\n%s", str(m1).c_str());
    auto det = get_determinant(m1);
    printf("Determinant = %f\n", det);

    m1.swap_to_row_major();
    printf("Matrix (swap to row-major)\n%s", str(m1).c_str());
    det = get_determinant(m1);
    printf("Determinant = %f\n", det);
}

void test_multiply()
{
    matrix_m<double> m1({{1.0, 2.0}, {2.0, -1.0}}, MAJOR::ROW);
    matrix_m<double> m2({{2.0, 0.0, -1.0}, {1.0, -1.0, 0.0}}, MAJOR::ROW);
    printf("Mat1\n%sMat2\n%s", str(m1).c_str(), str(m2).c_str());
    auto m1_times_m2 = m1 * m2;
    printf("Mat1*Mat2\n%s", str(m1_times_m2).c_str());
    matrix_m<double> m3({{4.0, -2.0, -1.0}, {3.0, 1.0, -2.0}}, MAJOR::ROW);
    assert(fequal_array(6, m1_times_m2.c, m3.c));

    // same result should be obtained with column-major
    m1.swap_to_col_major();
    m2.swap_to_col_major();
    m1_times_m2 = m1 * m2;
    m3.swap_to_col_major();
    assert(fequal_array(6, m1_times_m2.c, m3.c));
}

int main (int argc, char *argv[])
{
    test_constuctor();
    test_multiply();
    return 0;
}
