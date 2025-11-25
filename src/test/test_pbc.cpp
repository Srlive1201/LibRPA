#include "../core/pbc.h"
#include <cassert>

using namespace librpa_int;

static void test_is_gamma_point()
{
    assert(is_gamma_point(Vector3_Order<double>{0.0, 0.0, 0.0}));
    assert(!is_gamma_point(Vector3_Order<double>{0.3, 0.2, 0.1}));
    assert(is_gamma_point(Vector3_Order<int>{0, 0, 0}));
}

static void test_get_R_index()
{
    Vector3_Order<int> period{2, 2, 2};
    std::vector<Vector3_Order<int>> sc222 = construct_R_grid(period);
    for ( int i = 0; i != sc222.size(); i++ )
        printf("i=%d: x %d, y %d, z %d\n", i, sc222[i].x, sc222[i].y, sc222[i].z);
    Vector3_Order<int> R {-1, -1, -1};
    printf(" R: x %d, y %d, z %d\n", R.x, R.y, R.z);
    auto mR = -R;
    printf("-R: x %d, y %d, z %d\n", mR.x, mR.y, mR.z);
    assert(get_R_index(sc222, R) == 0);
    assert(get_R_index(sc222, Vector3_Order<int>{-1, 0, -1}) == 2);
    assert(get_R_index(sc222, Vector3_Order<int>{3, 3, -1}) < 0);
    assert(get_R_index(sc222, Vector3_Order<int>{3, 3, -1} % period) == 0);
}

static void test_periodic_boundary_data()
{
    PeriodicBoundaryData pbc;
    pbc.set_latvec_and_G({1, 2, 3, 4, 5, 6, 7, 8, 9},
                         {1, 2, 3, 4, 5, 6, 7, 8, 9});
}

int main (int argc, char *argv[])
{
    test_is_gamma_point();
    test_get_R_index();
    test_periodic_boundary_data();
    return 0;
}
