/*!
 * @file coulmat.h
 * @brief utilities for handling with Coulomb matrix
 */
#pragma once
#include "ri.h"

/*!
 * @brief perform Fourier transform of q/k-space Coulomb matrix to R-space
 * @param coulmat  q/k-space Coulomb matrix by atom-pair index
 * @param Rlist  vectors of R where the R-space Coulomb matrix should be computed
 * @param return_ordered_atom_pair  ensure that the returned type
 *        will always have keys of ordered atom pairs,
 *        i.e. there always exists map[I][J] and map[J][I] for I not equal to J
 * @return  atpair_R_mat_t
 */
atpair_R_mat_t FT_Vq(const atpair_k_cplx_mat_t &coulmat,
                     const int &n_k_points,
                     const vector<Vector3_Order<int>> &Rlist,
                     bool return_ordered_atom_pair);
