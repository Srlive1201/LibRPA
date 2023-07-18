/*!
 @file pbc.h
 @brief Utilities to deal with periodic boundary conditions
 */
#pragma once
#include <array>
#include <map>
#include <vector>
#include "vector3_order.h"
#include "matrix3.h"

// TODO: make it into a template
std::vector<Vector3_Order<int>> construct_R_grid(const Vector3_Order<int> &period);

//! Get the index of R in an Rlist. If R is not found in the list, return a negative number
int get_R_index(const std::vector<Vector3_Order<int>> &Rlist, const Vector3_Order<int> &R);

bool is_gamma_point(const Vector3_Order<double> &kpt);
bool is_gamma_point(const Vector3_Order<int> &kpt);

extern int kv_nmp[3];
//! lattice vectors as a 3D-matrix, each row as a lattice vector
extern Matrix3 latvec;
//! same as latvec, but a nested array for LibRI call
extern std::array<std::array<double, 3>, 3> lat_array;
//! reciprocal lattice vectors as a 3D-matrix, each row as a reciprocal vector
extern Matrix3 G;
extern std::vector<Vector3_Order<double>> klist;
extern std::vector<Vector3_Order<double>> kfrac_list;
extern std::vector<int> irk_point_id_mapping;
extern map<Vector3_Order<double>, vector<Vector3_Order<double>>> map_irk_ks;
extern Vector3<double> *kvec_c;
