#ifdef NDEBUG
#undef NDEBUG
#endif
#include "../qsgw/vxc_projection.h"
#include <cassert>
#include <cmath>
using namespace librpa_int;
using namespace librpa_int::qsgw;
void assert_close(cplxdb actual, cplxdb expected)
{
    assert(std::abs(actual - expected) < 1.e-14);
}
void test_abacus_nao_vxc_is_projected_to_the_fixed_state_basis()
{
    MeanField reference(1, 1, 2, 2, 1);
    const double inv_sqrt_two = 1.0 / std::sqrt(2.0);
    ComplexMatrix wfc(2, 2);
    wfc(0, 0) = inv_sqrt_two;
    wfc(0, 1) = cplxdb(0.0, inv_sqrt_two);
    wfc(1, 0) = cplxdb(0.0, inv_sqrt_two);
    wfc(1, 1) = inv_sqrt_two;
    reference.get_eigenvectors()[0][0][0] = wfc;

    Matz vxc_nao(2, 2);
    vxc_nao(0, 0) = 2.0;
    vxc_nao(0, 1) = cplxdb(1.0, 0.5);
    vxc_nao(1, 0) = cplxdb(1.0, -0.5);
    vxc_nao(1, 1) = 4.0;

    Matz coefficients(2, 2);
    for (int row = 0; row < 2; ++row)
    {
        for (int column = 0; column < 2; ++column)
        {
            coefficients(row, column) = wfc(row, column);
        }
    }
    const Matz expected =
        conj(coefficients) * vxc_nao * transpose(coefficients, false);
    const Matz actual =
        project_vxc_nao_to_fixed_basis(vxc_nao, reference, 0, 0);
    const Matz selected = prepare_vxc_in_fixed_state_basis(
        vxc_nao, VxcBasis::Nao, reference, 0, 0);

    for (int row = 0; row < 2; ++row)
    {
        for (int column = 0; column < 2; ++column)
        {
            assert_close(actual(row, column), expected(row, column));
            assert_close(selected(row, column), expected(row, column));
        }
    }
    assert(std::abs(actual(0, 1) - vxc_nao(0, 1)) > 1.0e-6);

    Matz vxc_state(2, 2);
    vxc_state(0, 0) = -1.0;
    vxc_state(0, 1) = cplxdb(0.2, 0.1);
    vxc_state(1, 0) = cplxdb(0.2, -0.1);
    vxc_state(1, 1) = 0.5;
    const Matz state_selected = prepare_vxc_in_fixed_state_basis(
        vxc_state, VxcBasis::State, reference, 0, 0);
    for (int row = 0; row < 2; ++row)
    {
        for (int column = 0; column < 2; ++column)
        {
            assert_close(state_selected(row, column), vxc_state(row, column));
        }
    }
}

int main()
{
    test_abacus_nao_vxc_is_projected_to_the_fixed_state_basis();
    return 0;
}
