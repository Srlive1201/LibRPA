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
    bool _lattice_set = false;
    bool _kgrids_kvec_set = false;
public:
    //! Lattice vectors as a 3D-matrix, each row as a lattice vector. Unit: Bohr
    Matrix3 latvec;
    //! Reciprocal lattice vectors as a 3D-matrix, each row as a reciprocal vector. Unit: 2pi/Bohr
    Matrix3 G;
    //! Lattice vectors as a nested array for LibRI call
    std::array<std::array<double, 3>, 3> latvec_array;

    //! Period of each lattice vector in the BvK cell
    Vector3_Order<int> period;

    std::vector<Vector3_Order<int>> Rlist;
    //! Full K-points in Cartesian coordinates
    std::vector<Vector3_Order<double>> klist;
    //! Full K-points in fractional coordinates, in units of G
    std::vector<Vector3_Order<double>> kfrac_list;

    //! Irreducible k-points in Cartesian coordinates
    std::vector<Vector3_Order<double>> klist_ibz;
    std::vector<double> kweight_ibz;
    std::vector<int> irk_point_id_mapping;
    //! Mapping of irreducible k-points to the list of k-points connected by symmetry
    std::map<Vector3_Order<double>, std::vector<Vector3_Order<double>>> map_irk_ks;

    PeriodicBoundaryData() = default;

    // Setting
    void set_latvec(const std::vector<double> &lat_mat);
    // Also set reciprocal lattice vectors by hand. The parsed rec_mat should be in Bohr^-1 unit.
    // No consistency check.
    void set_latvec_and_G(const std::vector<double> &lat_mat,
                          const std::vector<double> &rec_mat);
    void set_kgrids_kvec(int nk1, int nk2, int nk3, std::vector<double> &kvecs);
    bool is_set_finished();

    // Getting
    int get_R_index(const Vector3_Order<int> &R);
    int get_k_index_full(const Vector3_Order<double> &k);
};

// TODO: make it into a template
std::vector<Vector3_Order<int>> construct_R_grid(const Vector3_Order<int> &period);

//! Get the index of R in an Rlist. If R is not found in the list, return a negative number
int get_R_index(const std::vector<Vector3_Order<int>> &Rlist, const Vector3_Order<int> &R);

bool is_gamma_point(const Vector3_Order<double> &kpt);
bool is_gamma_point(const Vector3_Order<int> &kpt);

extern int kv_nmp[3];
//! lattice vectors as a 3D-matrix, each row as a lattice vector. Unit: Bohr
extern Matrix3 latvec;
//! same as latvec, but a nested array for LibRI call
extern std::array<std::array<double, 3>, 3> lat_array;
//! reciprocal lattice vectors as a 3D-matrix, each row as a reciprocal vector. Unit: 2pi/Bohr
extern Matrix3 G;
extern std::vector<Vector3_Order<double>> klist;
extern std::vector<Vector3_Order<double>> klist_ibz;
extern std::vector<Vector3_Order<double>> kfrac_list;
extern std::vector<int> irk_point_id_mapping;
extern std::map<Vector3_Order<double>, std::vector<Vector3_Order<double>>> map_irk_ks;
extern Vector3<double> *kvec_c;

int get_k_index_full(const Vector3_Order<double> &k);

}
