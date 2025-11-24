#include "pbc.h"
#include <algorithm>

namespace librpa_int {

std::vector<Vector3_Order<int>> construct_R_grid(const Vector3_Order<int> &period)
{
    // cout<<" begin to construct_R_grid"<<endl;
    std::vector<Vector3_Order<int>> R_grid;
    R_grid.clear();

    for (int x = -(period.x) / 2; x <= (period.x - 1) / 2; ++x)
        for (int y = -(period.y) / 2; y <= (period.y - 1) / 2; ++y)
            for (int z = -(period.z) / 2; z <= (period.z - 1) / 2; ++z)
                R_grid.push_back({x, y, z});

    return R_grid;
}

int get_R_index(const std::vector<Vector3_Order<int>> &Rlist, const Vector3_Order<int> &R)
{
    // v1: manual search
    // for ( int iR = 0; iR != Rlist.size(); iR++ )
    // {
    //     if (Rlist[iR] == R) return iR;
    // }
    // v2: use algorithm
    auto itr = std::find(Rlist.cbegin(), Rlist.cend(), R);
    if ( itr != Rlist.cend()) return distance(Rlist.cbegin(), itr);
    return -1;
}

bool is_gamma_point(const Vector3_Order<double> &kpt)
{
    double thres = 1.0e-5;
    return -thres < kpt.x && kpt.x < thres
        && -thres < kpt.y && kpt.y < thres
        && -thres < kpt.z && kpt.z < thres;
}

bool is_gamma_point(const Vector3_Order<int> &kpt_int)
{
    return kpt_int.x == 0 && kpt_int.y == 0 && kpt_int.z == 0;
}

int kv_nmp[3] = {1, 1, 1};
Vector3<double> *kvec_c;
std::vector<Vector3_Order<double>> klist;
std::vector<Vector3_Order<double>> klist_ibz;
std::vector<Vector3_Order<double>> kfrac_list;
std::vector<int> irk_point_id_mapping;
std::map<Vector3_Order<double>, std::vector<Vector3_Order<double>>> map_irk_ks;
Matrix3 latvec;
std::array<std::array<double, 3>, 3> lat_array;
Matrix3 G;

int get_k_index_full(const Vector3_Order<double> &k)
{
    auto it = std::find(klist.cbegin(), klist.cend(), k);
    if ( it != klist.cend()) return distance(klist.cbegin(), it);
    return -1;
}

}
