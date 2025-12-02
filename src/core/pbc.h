/*!
 @file pbc.h
 @brief Utilities to deal with periodic boundary conditions
 */
#pragma once
#include <array>
#include <map>
#include <vector>

#include "../math/matrix3.h"
#include "../math/vector3_order.h"

namespace librpa_int {

class PeriodicBoundaryData
{
private:
    bool lattice_reset_ = false;
    bool kgrid_no_symmetry_ = true;

public:
    //! Lattice vectors as a 3D-matrix, each row as a lattice vector. Unit: Bohr
    Matrix3 latvec;
    //! Reciprocal lattice vectors as a 3D-matrix, each row as a reciprocal vector. Unit: 2pi/Bohr
    Matrix3 G;
    //! Lattice vectors as a nested array for LibRI call
    std::array<std::array<double, 3>, 3> latvec_array;

    //! Period of each lattice vector in the BvK cell
    Vector3_Order<int> period;
    //! Period as an array for LibRI call
    std::array<int, 3> period_array;

    std::vector<Vector3_Order<int>> Rlist;
    //! Full K-points in Cartesian coordinates
    std::vector<Vector3_Order<double>> klist;
    //! Full K-points in fractional coordinates
    std::vector<Vector3_Order<double>> kfrac_list;

    //! Irreducible k-points in Cartesian coordinates
    std::vector<Vector3_Order<double>> klist_ibz;
    //! Scalar weight of each irreducible k-points
    std::vector<double> kweight_ibz;
    //! Same as kweight_ibz, but indexed directly by the k-point
    std::map<Vector3_Order<double>, double> map_ibzk_weight;
    std::vector<int> irk_point_id_mapping;
    std::vector<int> isymops;
    //! Mapping of irreducible k-points to the list of k-points connected by symmetry
    std::map<Vector3_Order<double>, std::vector<Vector3_Order<double>>> map_irk_ks;

    PeriodicBoundaryData();

    //! Setting lattice vectors. latt_mat should be in Bohr unit.
    void set_latvec(const std::vector<double> &lat_mat);
    //! Set lattice and reciprocal lattice vectors. The parsed rec_mat should be in Bohr^-1 unit. No consistency check.
    void set_latvec_and_G(const std::vector<double> &latt_mat,
                          const std::vector<double> &recp_mat);
    //! Set the full k-points list. kvecs should be in Bohr^-1 unit.
    void set_kgrids_kvec(int nk1, int nk2, int nk3, std::vector<double> &kvecs);
    //! Set the mapping of full k-points to irreducible k-points
    void set_ibz_mapping(const std::vector<int> &irk_point_id_mapping_in,
                         const std::vector<int> &isymops_in = {});

    // Getting
    int get_R_index(const Vector3_Order<int> &R) const;
    int get_k_index_full(const Vector3_Order<double> &k) const;
    int get_n_cells_bvk() const { return period.x * period.y * period.z; }
    bool is_latt_set() const { return lattice_reset_; }
};

// TODO: make it into a template
std::vector<Vector3_Order<int>> construct_R_grid(const Vector3_Order<int> &period);

//! Get the index of R in an Rlist. If R is not found in the list, return a negative number
int get_R_index(const std::vector<Vector3_Order<int>> &Rlist, const Vector3_Order<int> &R);

bool is_gamma_point(const Vector3_Order<double> &kpt, double thres = 1.0e-5);
bool is_gamma_point(const Vector3_Order<int> &kpt);

// extern int kv_nmp[3];
// //! lattice vectors as a 3D-matrix, each row as a lattice vector. Unit: Bohr
// extern Matrix3 latvec;
// //! same as latvec, but a nested array for LibRI call
// extern std::array<std::array<double, 3>, 3> lat_array;
// //! reciprocal lattice vectors as a 3D-matrix, each row as a reciprocal vector. Unit: 2pi/Bohr
// extern Matrix3 G;
// extern std::vector<Vector3_Order<double>> klist;
// extern std::vector<Vector3_Order<double>> klist_ibz;
// extern std::vector<Vector3_Order<double>> kfrac_list;
// extern std::vector<int> irk_point_id_mapping;
// extern std::map<Vector3_Order<double>, std::vector<Vector3_Order<double>>> map_irk_ks;
// extern Vector3<double> *kvec_c;

int get_k_index_full(const Vector3_Order<double> &k);

}
