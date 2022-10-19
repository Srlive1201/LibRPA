/*
 * @file coulmat.h
 * @brief utilities for handling with Coulomb matrix
 */
#pragma once
#include "ri.h"

/*
 * @brief perform Fourier transform of q/k-space Coulomb matrix to R-space
 * @param coulmat  q/k-space Coulomb matrix by atom-pair index
 * @param Rlist  vectors of R where the R-space Coulomb matrix should be computed
 * @return  atpair_R_mat_t
 */
atpair_R_mat_t FT_Vq(const atpair_k_cplx_mat_t &coulmat, vector<Vector3_Order<int>> Rlist);

