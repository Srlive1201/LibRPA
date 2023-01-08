#include "../src/pbc.h"
#include <cassert>

int main (int argc, char *argv[])
{
    Vector3_Order<int> period{2, 2, 2};
    vector<Vector3_Order<int>> sc222 = construct_R_grid(period);
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
    return 0;
}
