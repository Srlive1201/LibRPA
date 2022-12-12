/*!
 @file pbc.h
 @brief Utilities to deal with periodic boundary conditions
 */
#pragma once
#include <vector>
#include "vector3_order.h"
using std::vector;

// TODO: make it into a template
vector<Vector3_Order<int>> construct_R_grid(const Vector3_Order<int> & period);

//! Get the index of R in an Rlist. If R is not found in the list, return a negative number
int get_R_index(vector<Vector3_Order<int>> Rlist, const Vector3_Order<int> &R);

bool is_gamma_point(const Vector3_Order<double> &kpt);
