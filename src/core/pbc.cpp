#include "pbc.h"

#include <algorithm>
#include <numeric>
#include <stdexcept>

#include "../utils/constants.h"

namespace librpa_int {

void PeriodicBoundaryData::set_latvec(const std::vector<double> &lat_mat)
{
    if (lat_mat.size() != 9)
        throw std::runtime_error("Invalid lattice vectors input");

    this->latvec.e11 = lat_mat[0];
    this->latvec.e12 = lat_mat[1];
    this->latvec.e13 = lat_mat[2];
    this->latvec.e21 = lat_mat[3];
    this->latvec.e22 = lat_mat[4];
    this->latvec.e23 = lat_mat[5];
    this->latvec.e31 = lat_mat[6];
    this->latvec.e32 = lat_mat[7];
    this->latvec.e33 = lat_mat[8];

    this->latvec_array[0] = {this->latvec.e11,this->latvec.e12,this->latvec.e13};
    this->latvec_array[1] = {this->latvec.e21,this->latvec.e22,this->latvec.e23};
    this->latvec_array[2] = {this->latvec.e31,this->latvec.e32,this->latvec.e33};

    this->G = this->latvec.Inverse().Transpose();

    this->_lattice_set = true;
}

void PeriodicBoundaryData::set_latvec_and_G(const std::vector<double> &lat_mat,
                                            const std::vector<double> &rec_mat)
{
    this->set_latvec(lat_mat);

    if (rec_mat.size() != 9)
        throw std::runtime_error("Invalid reciprocal lattice vectors input");

    this->G.e11 = rec_mat[0];
    this->G.e12 = rec_mat[1];
    this->G.e13 = rec_mat[2];
    this->G.e21 = rec_mat[3];
    this->G.e22 = rec_mat[4];
    this->G.e23 = rec_mat[5];
    this->G.e31 = rec_mat[6];
    this->G.e32 = rec_mat[7];
    this->G.e33 = rec_mat[8];
    // Scale to 2pi/bohr unit
    this->G /= TWO_PI;

    this->_lattice_set = true;
}

void PeriodicBoundaryData::set_kgrids_kvec(int nk1, int nk2, int nk3, std::vector<double> &kvecs)
{
    this->period = {nk1, nk2, nk3};

    this->Rlist = construct_R_grid(this->period);
    klist.clear();
    kfrac_list.clear();
    int n_k_points = nk1 * nk2 * nk3;
    for (int ik = 0; ik < n_k_points; ik++)
    {
        double kx = kvecs[ik * 3];
        double ky = kvecs[ik * 3 + 1];
        double kz = kvecs[ik * 3 + 2];
        Vector3_Order<double> kvec{kx, ky, kz};
        klist.emplace_back(kvec);
        kfrac_list.emplace_back(latvec * kvec);
    }

    // Initialize the irreducible k points: same as the full ones
    klist_ibz.clear();
    klist_ibz.insert(klist_ibz.begin(), klist.begin(), klist.end());

    kweight_ibz.clear();
    kweight_ibz.resize(n_k_points, 1.0 / n_k_points);

    irk_point_id_mapping.clear();
    std::iota(irk_point_id_mapping.begin(), irk_point_id_mapping.end(), 0);

    map_irk_ks.clear();
    for (const auto &k: klist)
    {
        map_irk_ks[k].emplace_back(k);
    }

    this->_kgrids_kvec_set = true;
}

inline bool PeriodicBoundaryData::is_set_finished()
{
    return this->_kgrids_kvec_set && this->_lattice_set;
}

inline int PeriodicBoundaryData::get_R_index(const Vector3_Order<int> &R)
{
    return librpa_int::get_R_index(this->Rlist, R);
}

int PeriodicBoundaryData::get_k_index_full(const Vector3_Order<double> &k)
{
    auto it = std::find(this->klist.cbegin(), this->klist.cend(), k);
    if ( it != this->klist.cend()) return distance(this->klist.cbegin(), it);
    return -1;
}

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
