#include "../math/matrix_m.h"
#include "testutils.h"

using namespace std;
using namespace librpa_int;

static void test_constuctor()
{
    Matd m1(2, 3, MAJOR::COL);
    m1.zero_out();
    m1(0, 1) = 1.0;
    m1(1, 0) = 2.0;
    printf("Matrix (set column-major)\n%s", str(m1).c_str());
    auto det = get_determinant(m1);
    printf("Determinant = %f\n", det);
    assert(fequal(det, -2.0));

    m1.swap_to_row_major();
    printf("Matrix (swap to row-major)\n%s", str(m1).c_str());
    det = get_determinant(m1);
    printf("Determinant = %f\n", det);
    assert(fequal(det, -2.0));
}

static void test_minmax()
{
    int ir = 0, ic = 0;
    Matd m1({{1.0, 2.0}, {3.0, -4.0}}, MAJOR::ROW);
    assert(fequal(m1.max(), 3.0));
    assert(fequal(m1.min(), -4.0));
    m1.argmin(ir, ic);
    assert( ir == 1 && ic == 1);
    m1.argmax(ir, ic);
    assert( ir == 1 && ic == 0);

    // complex matrix
    Matz mc1({{1.0, 2.0}, {3.0, -4.0}}, MAJOR::ROW);
    assert(fequal(mc1.max(), 4.0));
    assert(fequal(mc1.min(), 1.0));
    mc1.argmin(ir, ic);
    assert( ir == 0 && ic == 0);
    mc1.argmax(ir, ic);
    assert( ir == 1 && ic == 1);
}

static void test_power_hemat()
{
    cout << "Run test_power_hemat" << endl;
    const int n = 8;
    auto hemat = random_he_selected_ev<double>(n, {1.0, 2.25, 4.0, 9.0}, MAJOR::ROW);
    cout << "hemat" << endl << hemat;
    const auto hemat_back = hemat.copy();
    auto sqrt_hemat = power_hemat(hemat, 0.5, false, false);
    cout << "hemat after power_hemat" << endl << hemat;
    assert(hemat == hemat_back);
    assert(hemat == sqrt_hemat * sqrt_hemat);
    hemat.swap_to_col_major();
    sqrt_hemat = power_hemat(hemat, 0.5, false, false);
    assert(hemat == sqrt_hemat * sqrt_hemat);
}

static void test_multiply()
{
    Matd m1({{1.0, 2.0}, {2.0, -1.0}}, MAJOR::ROW);
    Matd m2({{2.0, 0.0, -1.0}, {1.0, -1.0, 0.0}}, MAJOR::ROW);
    printf("Mat1\n%sMat2\n%s", str(m1).c_str(), str(m2).c_str());
    auto m1_times_m2 = m1 * m2;
    printf("Mat1*Mat2\n%s", str(m1_times_m2).c_str());
    Matd m3({{4.0, -2.0, -1.0}, {3.0, 1.0, -2.0}}, MAJOR::ROW);
    assert(fequal_array(6, m1_times_m2.ptr(), m3.ptr()));

    // same result should be obtained with column-major
    m1.swap_to_col_major();
    m2.swap_to_col_major();
    m1_times_m2 = m1 * m2;
    m3.swap_to_col_major();
    assert(fequal_array(6, m1_times_m2.ptr(), m3.ptr()));
}

static void test_conjugate()
{
    Matz m({{{1.0, -2.0}, {2.0, -1.0}}, {{2.0, 1.0}, -1.0}});
    Matz m_conj_ref({{{1.0, 2.0}, {2.0, 1.0}}, {{2.0, -1.0}, -1.0}});
    m.conj();
    assert(m == m_conj_ref);
}

int main (int argc, char *argv[])
{
    test_constuctor();
    test_multiply();
    test_minmax();
    test_power_hemat();
    test_conjugate();
    return 0;
}
