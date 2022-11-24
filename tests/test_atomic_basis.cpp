#include <cassert>
#include <unordered_map>
#include <utility>
#include "../src/atomic_basis.h"
#include "testutils.h"

void test_constuctor()
{
    LIBRPA::AtomicBasis ab;
    assert(ab.n_atoms == 0);
    assert(ab.nb_total == 0);

    ab.set({5, 4, 3});
    assert(ab.n_atoms == 3);
    assert(ab.nb_total == 12);

    LIBRPA::AtomicBasis ab_from_map({0, 1, 1}, {{0, 4}, {1, 8}});
    assert(ab_from_map.n_atoms == 3);
    assert(ab_from_map.nb_total == (4 + 2 * 8));
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

int main (int argc, char *argv[])
{
    test_constuctor();
    test_indexing();
    return 0;
}
