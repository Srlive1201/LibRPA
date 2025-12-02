#include "pbc.h"

#include <algorithm>
#include <numeric>
#include <stdexcept>

#include "../utils/constants.h"
#include "../utils/error.h"

namespace librpa_int {

static const double a = 1e5;  // 10000 nm
static const double inva = 1.0 / a;

// Initialized as a huge box to emulate isolated case
PeriodicBoundaryData::PeriodicBoundaryData():
    lattice_reset_(false),
    kgrid_no_symmetry_(true),
    latvec(a, 0, 0, 0, a, 0, 0, 0, a),
    G(inva, 0, 0, 0, inva, 0, 0, 0, inva),
    period(1, 1, 1),
    period_array({1, 1, 1}),
    Rlist({{0, 0, 0}}),
    klist({{0, 0, 0}}),
    kfrac_list({{0, 0, 0}})
{
    this->latvec_array[0] = {this->latvec.e11,this->latvec.e12,this->latvec.e13};
    this->latvec_array[1] = {this->latvec.e21,this->latvec.e22,this->latvec.e23};
    this->latvec_array[2] = {this->latvec.e31,this->latvec.e32,this->latvec.e33};

    // Gamma-only
    this->set_ibz_mapping({0}, {0});
}

void PeriodicBoundaryData::set_latvec(const std::vector<double> &latt_mat)
{
    if (latt_mat.size() != 9)
        throw std::runtime_error("Invalid lattice vectors input");

    this->latvec.e11 = latt_mat[0];
    this->latvec.e12 = latt_mat[1];
    this->latvec.e13 = latt_mat[2];
    this->latvec.e21 = latt_mat[3];
    this->latvec.e22 = latt_mat[4];
    this->latvec.e23 = latt_mat[5];
    this->latvec.e31 = latt_mat[6];
    this->latvec.e32 = latt_mat[7];
    this->latvec.e33 = latt_mat[8];

    this->latvec_array[0] = {this->latvec.e11,this->latvec.e12,this->latvec.e13};
    this->latvec_array[1] = {this->latvec.e21,this->latvec.e22,this->latvec.e23};
    this->latvec_array[2] = {this->latvec.e31,this->latvec.e32,this->latvec.e33};

    this->G = this->latvec.Inverse().Transpose();

    lattice_reset_ = true;
}

void PeriodicBoundaryData::set_latvec_and_G(const std::vector<double> &latt_mat,
                                            const std::vector<double> &recp_mat)
{
    this->set_latvec(latt_mat);

    if (recp_mat.size() != 9)
        throw std::runtime_error("Invalid reciprocal lattice vectors input");

    this->G.e11 = recp_mat[0];
    this->G.e12 = recp_mat[1];
    this->G.e13 = recp_mat[2];
    this->G.e21 = recp_mat[3];
    this->G.e22 = recp_mat[4];
    this->G.e23 = recp_mat[5];
    this->G.e31 = recp_mat[6];
    this->G.e32 = recp_mat[7];
    this->G.e33 = recp_mat[8];
    this->G /= TWO_PI;
}

void PeriodicBoundaryData::set_kgrids_kvec(int nk1, int nk2, int nk3, std::vector<double> &kvecs)
{
    this->period = {nk1, nk2, nk3};
    this->period_array = {nk1, nk2, nk3};

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
        kvec /= TWO_PI;
        klist.emplace_back(kvec);
        kfrac_list.emplace_back(latvec * kvec);
    }

    // Initialize the irreducible k points: same as the full ones
    std::vector<int> irk_point_id_mapping_in(n_k_points);
    std::iota(irk_point_id_mapping_in.begin(), irk_point_id_mapping_in.end(), 0);
    this->set_ibz_mapping(irk_point_id_mapping_in);
}

void PeriodicBoundaryData::set_ibz_mapping(const std::vector<int> &irk_point_id_mapping_in,
                                           const std::vector<int> &isymops_in)
{
    const size_t &n_k_points = this->klist.size();
    const double step = 1.0 / n_k_points;

    if (irk_point_id_mapping_in.size() != n_k_points)
    {
        throw LIBRPA_RUNTIME_ERROR("irk_point_id_mapping size not equal to n_k_points: " + std::to_string(n_k_points));
    }
    // TODO: implement using symmetry operation mapping
    this->isymops = isymops_in;

    irk_point_id_mapping = irk_point_id_mapping_in;

    klist_ibz.clear();
    map_ibzk_weight.clear();
    map_irk_ks.clear();
    for (size_t ik = 0; ik < n_k_points; ik++)
    {
        int ik_ibz = irk_point_id_mapping[ik];
        const auto &k = klist[ik];
        const auto &k_ibz = klist[ik_ibz];
        auto it = map_ibzk_weight.find(k_ibz);
        if (it == map_ibzk_weight.end())
        {
            klist_ibz.emplace_back(k_ibz);
            map_ibzk_weight.emplace(k_ibz, step);
        }
        else
        {
            it->second += step;
        }
        map_irk_ks[k_ibz].emplace_back(k);
    }

    kweight_ibz.clear();
    for (const auto &k_ibz: klist_ibz)
    {
        kweight_ibz.emplace_back(map_ibzk_weight.at(k_ibz));
    }

    kgrid_no_symmetry_ = klist_ibz.size() == klist.size();
}

int PeriodicBoundaryData::get_R_index(const Vector3_Order<int> &R) const
{
    return librpa_int::get_R_index(this->Rlist, R);
}

static int get_k_index_(const std::vector<Vector3_Order<double>> &klist, const Vector3_Order<double> &k)
{
    auto it = std::find(klist.cbegin(), klist.cend(), k);
    if (it != klist.cend()) return distance(klist.cbegin(), it);
    return -1;
}

int PeriodicBoundaryData::get_k_index_full(const Vector3_Order<double> &k) const
{
    return get_k_index_(this->klist, k);
}

int PeriodicBoundaryData::get_k_index_ibz(const Vector3_Order<double> &k) const
{
    return get_k_index_(this->klist_ibz, k);
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
    auto itr = std::find(Rlist.cbegin(), Rlist.cend(), R);
    if ( itr != Rlist.cend()) return distance(Rlist.cbegin(), itr);
    return -1;
}

bool is_gamma_point(const Vector3_Order<double> &kpt, double thres)
{
    return -thres < kpt.x && kpt.x < thres
        && -thres < kpt.y && kpt.y < thres
        && -thres < kpt.z && kpt.z < thres;
}

bool is_gamma_point(const Vector3_Order<int> &kpt_int)
{
    return kpt_int.x == 0 && kpt_int.y == 0 && kpt_int.z == 0;
}

// int kv_nmp[3] = {1, 1, 1};
// Vector3<double> *kvec_c;
// std::vector<Vector3_Order<double>> klist;
// std::vector<Vector3_Order<double>> klist_ibz;
// std::vector<Vector3_Order<double>> kfrac_list;
// std::vector<int> irk_point_id_mapping;
// std::map<Vector3_Order<double>, std::vector<Vector3_Order<double>>> map_irk_ks;
// Matrix3 latvec;
// std::array<std::array<double, 3>, 3> lat_array;
// Matrix3 G;

}
