/*!
 @file pbc.h
 @brief Utilies to deal with periodic boundary conditions
 */
#pragma once
#include <vector>
#include "vector3_order.h"

// TODO: make it into a template
vector<Vector3_Order<int>> construct_R_grid(const Vector3_Order<int> & period);
