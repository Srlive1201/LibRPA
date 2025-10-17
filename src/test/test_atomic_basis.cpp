#include "../atomic_basis.h"

#include "testutils.h"

#include <cassert>
#include <unordered_map>
#include <utility>

void test_constuctor()
{
    LIBRPA::AtomicBasis ab;
    assert(ab.n_atoms == 0);
    assert(ab.nb_total == 0);

    ab.set({5, 4, 3});
    assert(ab.n_atoms == 3);
    assert(ab.nb_total == 12);

    // unordered iatom-nbasis set
    // 3 atoms 2, 5, 4
    ab.set({{1, 5}, {0, 2}, {4, 4}});
    assert(ab.n_atoms == 3);
    assert(ab.nb_total == 11);
    assert(8 == ab.get_global_index(2, 1));
    assert(4 == ab.get_global_index(1, 2));

    LIBRPA::AtomicBasis ab_from_doublemap({0, 1, 1}, {{0, 4}, {1, 8}});
    assert(ab_from_doublemap.n_atoms == 3);
    assert(ab_from_doublemap.nb_total == (4 + 2 * 8));
}

void test_indexing()
{
    LIBRPA::AtomicBasis ab({4, 4, 8, 8, 5});
    assert(ab.nb_total == 29);
    std::unordered_map<int, std::pair<int, int>> map_go_I_lo
        {
            {0, {0, 0}},
            {4, {1, 0}},
            {6, {1, 2}},
            {14, {2, 6}},
            {24, {4, 0}},
            {28, {4, 4}},
        };
    int I, ilo, igo;
    for (const auto& go_I_lo: map_go_I_lo)
    {
        igo = go_I_lo.first;
        ab.get_local_index(igo, I, ilo);
        assert(I == go_I_lo.second.first);
        assert(ilo == go_I_lo.second.second);
        assert(igo == ab.get_global_index(I, ilo));
    }
}

void test_get_2d_indices()
{
    LIBRPA::AtomicBasis ab({1, 2, 3});
    // (0, 1), (0, 2), (1, 0), (2, 0)
    const auto id_col_fast = LIBRPA::get_2d_mat_indices_atpair(ab, ab, {{0, 1}, {1, 0}}, false);
    assert(id_col_fast.size() == 4);
    assert(equal_pair(id_col_fast[0], {0, 1}));
    assert(equal_pair(id_col_fast[1], {0, 2}));
    assert(equal_pair(id_col_fast[2], {1, 0}));
    assert(equal_pair(id_col_fast[3], {2, 0}));
    // (1, 3), (2, 3), (1, 4), (2, 4), (1, 5), (2, 5)
    const auto id_row_fast = LIBRPA::get_2d_mat_indices_atpair(ab, ab, {{1, 2}}, true);
    assert(id_row_fast.size() == 6);
    assert(equal_pair(id_row_fast[0], {1, 3}));
    assert(equal_pair(id_row_fast[1], {2, 3}));
    assert(equal_pair(id_row_fast[2], {1, 4}));
    assert(equal_pair(id_row_fast[3], {2, 4}));
    assert(equal_pair(id_row_fast[4], {1, 5}));
    assert(equal_pair(id_row_fast[5], {2, 5}));
}

void test_get_1d_indices()
{
    LIBRPA::AtomicBasis ab({1, 2, 3});
    // 2D: (1, 3), (2, 3), (1, 4), (2, 4), (1, 5), (2, 5)
    // row-major: 9, 15, 10, 16, 11, 17
    // col-major: 19, 20, 25, 26, 31, 32
    const auto id_rfast_rmajor = LIBRPA::get_1d_mat_indices_atpair(ab, ab, {{1, 2}}, true, true);
    const std::vector<std::size_t> ref_rfast_rmajor({9, 15, 10, 16, 11, 17});
    assert(id_rfast_rmajor.size() == 6);
    assert(equal_array(6, id_rfast_rmajor.data(), ref_rfast_rmajor.data()));

    const auto id_rfast_cmajor = LIBRPA::get_1d_mat_indices_atpair(ab, ab, {{1, 2}}, true, false);
    const std::vector<std::size_t> ref_rfast_cmajor({19, 20, 25, 26, 31, 32});
    assert(id_rfast_cmajor.size() == 6);
    assert(equal_array(6, id_rfast_cmajor.data(), ref_rfast_cmajor.data()));
}

int main (int argc, char *argv[])
{
    test_constuctor();
    test_indexing();
    test_get_2d_indices();
    test_get_1d_indices();
    return 0;
}
