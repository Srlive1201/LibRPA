#include <cassert>
#include <cmath>
#include <complex>
#include <map>

#include "../core/epsilon.h"

namespace
{

using namespace librpa_int;

void assert_close(const std::complex<double> actual, const std::complex<double> expected,
                  const double tolerance = 1e-12)
{
    assert(std::abs(actual - expected) <= tolerance);
}

void test_dense_symmetry_collection_hermitizes_diagonal_atom_blocks()
{
    symmetry_atom_block_matrix_map_t blocks;

    ComplexMatrix block_00(2, 2);
    block_00(0, 0) = {1.0, 0.4};
    block_00(0, 1) = {2.0, 3.0};
    block_00(1, 0) = {4.0, -1.0};
    block_00(1, 1) = {5.0, -0.2};
    blocks[0][0] = block_00;

    ComplexMatrix block_01(2, 1);
    block_01(0, 0) = {6.0, 7.0};
    block_01(1, 0) = {8.0, 9.0};
    blocks[0][1] = block_01;

    ComplexMatrix block_11(1, 1);
    block_11(0, 0) = {10.0, 2.0};
    blocks[1][1] = block_11;

    const std::map<atom_t, size_t> atom_nabf{{0, 2}, {1, 1}};
    const auto dense =
        build_dense_symmetry_hermitian_matrix_from_local_blocks(blocks, atom_nabf, true);

    assert(dense.nr == 3);
    assert(dense.nc == 3);
    assert_close(dense(0, 0), {1.0, 0.0});
    assert_close(dense(0, 1), {3.0, 2.0});
    assert_close(dense(1, 0), {3.0, -2.0});
    assert_close(dense(1, 1), {5.0, 0.0});

    assert_close(dense(0, 2), {6.0, 7.0});
    assert_close(dense(1, 2), {8.0, 9.0});
    assert_close(dense(2, 0), {6.0, -7.0});
    assert_close(dense(2, 1), {8.0, -9.0});
    assert_close(dense(2, 2), {10.0, 0.0});
}

void test_dense_symmetry_collection_preserves_generic_q_diagonal_blocks()
{
    symmetry_atom_block_matrix_map_t blocks;
    ComplexMatrix block(1, 1);
    block(0, 0) = {2.1, -0.7};
    blocks[0][0] = block;

    const std::map<atom_t, size_t> atom_nabf{{0, 1}};
    const auto dense =
        build_dense_symmetry_hermitian_matrix_from_local_blocks(blocks, atom_nabf, false);

    assert_close(dense(0, 0), {2.1, -0.7});
}

}  // namespace

int main()
{
    test_dense_symmetry_collection_hermitizes_diagonal_atom_blocks();
    test_dense_symmetry_collection_preserves_generic_q_diagonal_blocks();
    return 0;
}
