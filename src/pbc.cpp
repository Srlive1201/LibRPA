#include "pbc.h"

vector<Vector3_Order<int>> construct_R_grid(const Vector3_Order<int> & period)
{
    // cout<<" begin to construct_R_grid"<<endl;
    vector<Vector3_Order<int>> R_grid;
    R_grid.clear();

    for (int x = -(period.x) / 2; x <= (period.x - 1) / 2; ++x)
        for (int y = -(period.y) / 2; y <= (period.y - 1) / 2; ++y)
            for (int z = -(period.z) / 2; z <= (period.z - 1) / 2; ++z)
                R_grid.push_back({x, y, z});

    return R_grid;
}
