#include "pbc.h"
#include <algorithm>

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

int get_R_index(vector<Vector3_Order<int>> Rlist, const Vector3_Order<int> &R)
{
    // v1: manual search
    // for ( int iR = 0; iR != Rlist.size(); iR++ )
    // {
    //     if (Rlist[iR] == R) return iR;
    // }
    // v2: use algorithm
    vector<Vector3_Order<int>>::iterator itr = std::find(Rlist.begin(), Rlist.end(), R);
    if ( itr != Rlist.end()) return distance(Rlist.begin(), itr);
    return -1;
}

bool is_gamma_point(const Vector3_Order<double> &kpt)
{
    double thres = 1.0e-5;
    return (-thres < kpt.x < thres) && (-thres < kpt.y < thres) && (-thres < kpt.z < thres);
}

bool is_gamma_point(const Vector3_Order<int> &kpt_int)
{
    return kpt_int.x == 0 && kpt_int.y == 0 && kpt_int.z == 0;
}
