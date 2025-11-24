#include "../core/atomic_basis.h"

// #include "../utils/stl_io_helper.h"
#include "testutils.h"

#include <cassert>
#include <unordered_map>
#include <utility>

using namespace librpa_int;

void test_constuctor()
{
    librpa_int::AtomicBasis ab;
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

    librpa_int::AtomicBasis ab_from_doublemap({0, 1, 1}, {{0, 4}, {1, 8}});
    assert(ab_from_doublemap.n_atoms == 3);
    assert(ab_from_doublemap.nb_total == (4 + 2 * 8));
}

void test_indexing()
{
    librpa_int::AtomicBasis ab({4, 4, 8, 8, 5});
    assert(ab.nb_total == 29);
    std::unordered_map<size_t, std::pair<int, int>> map_go_I_lo
        {
            {0, {0, 0}},
            {4, {1, 0}},
            {6, {1, 2}},
            {14, {2, 6}},
            {24, {4, 0}},
            {28, {4, 4}},
        };
    int I, ilo;
    size_t igo;
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
    librpa_int::AtomicBasis ab({1, 2, 3});
    {
        // column fast case
        const auto id = librpa_int::get_2d_mat_indices_atpair(ab, ab, {{0, 1}, {1, 0}}, false, false);
        assert(id.size() == 4);
        const std::vector<atpair_t> ref = {{0, 1}, {0, 2}, {1, 0}, {2, 0}};
        for (int i = 0; i < 4; i++) assert(equal_pair(id[i], ref[i]));
    }
    {
        // row fast case
        const auto id = librpa_int::get_2d_mat_indices_atpair(ab, ab, {{1, 2}}, true, false);
        assert(id.size() == 6);
        const std::vector<atpair_t> ref = {{1, 3}, {2, 3}, {1, 4}, {2, 4}, {1, 5}, {2, 5}};
        for (int i = 0; i < 6; i++) assert(equal_pair(id[i], ref[i]));
    }
    {
        // another row fast case
        const std::vector<atpair_t> ref_nosort = {
            {1, 1}, {2, 1}, {1, 2}, {2, 2}, // (1, 1)
            {3, 1}, {4, 1}, {5, 1}, {3, 2}, {4, 2}, {5, 2}, // (2, 1)
        };
        const std::vector<atpair_t> ref_sort = {
            {1, 1}, {2, 1}, {3, 1}, {4, 1}, {5, 1}, // second column
            {1, 2}, {2, 2}, {3, 2}, {4, 2}, {5, 2}, // third column
        };
        // without sort
        auto id = librpa_int::get_2d_mat_indices_atpair(ab, ab, {{1, 1}, {2, 1}}, true, false);
        assert(id.size() == 10);
        for (int i = 0; i < 10; i++) assert(equal_pair(id[i], ref_nosort[i]));
        // sort
        id = librpa_int::get_2d_mat_indices_atpair(ab, ab, {{1, 1}, {2, 1}}, true, true);
        assert(id.size() == 10);
        for (int i = 0; i < 10; i++) assert(equal_pair(id[i], ref_sort[i]));
    }
    {
        // almost same as above, but column fast
        const std::vector<atpair_t> ref_nosort = {
            {1, 1}, {1, 2}, {2, 1}, {2, 2}, // (1, 1)
            {3, 1}, {3, 2}, {4, 1}, {4, 2}, {5, 1}, {5, 2}, // (2, 1)
        };
        const std::vector<atpair_t> ref_sort = {
            {1, 1}, {1, 2}, {2, 1}, {2, 2}, {3, 1}, {3, 2}, {4, 1}, {4, 2}, {5, 1}, {5, 2},
        };
        // without sort
        auto id = librpa_int::get_2d_mat_indices_atpair(ab, ab, {{1, 1}, {2, 1}}, false, false);
        assert(id.size() == 10);
        for (int i = 0; i < 10; i++) assert(equal_pair(id[i], ref_nosort[i]));
        // sort
        id = librpa_int::get_2d_mat_indices_atpair(ab, ab, {{1, 1}, {2, 1}}, false, true);
        assert(id.size() == 10);
        for (int i = 0; i < 10; i++) assert(equal_pair(id[i], ref_sort[i]));
    }

}

void test_get_1d_indices()
{
    librpa_int::AtomicBasis ab({1, 2, 3});
    // 2D: (1, 3), (2, 3), (1, 4), (2, 4), (1, 5), (2, 5)
    // row-major: 9, 15, 10, 16, 11, 17
    // col-major: 19, 20, 25, 26, 31, 32
    const auto id_rfast_rmajor = librpa_int::get_1d_mat_indices_atpair(ab, ab, {{1, 2}}, true, true);
    const std::vector<std::size_t> ref_rfast_rmajor({9, 15, 10, 16, 11, 17});
    assert(id_rfast_rmajor.size() == 6);
    assert(equal_array(6, id_rfast_rmajor.data(), ref_rfast_rmajor.data()));

    const auto id_rfast_cmajor = librpa_int::get_1d_mat_indices_atpair(ab, ab, {{1, 2}}, true, false);
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
