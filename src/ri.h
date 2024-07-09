/*!
 @file ri.h
 @brief Utilies related to resolution of identity
 */
#ifndef RI_H
#define RI_H
#include <cstring>
#include <map>
#include <memory>
#include <utility>
#include "vector3_order.h"
#include "atoms.h"
#include "complexmatrix.h"
using std::map;

extern int n_irk_points;
extern int natom;
extern int ncell;
extern std::vector<Vector3_Order<double>> irk_points;
extern map<Vector3_Order<double>, double> irk_weight;
extern map<atom_t, size_t> atom_nw;
extern map<atom_t, size_t> atom_mu;
extern map<atom_t, size_t> atom_nw_loc;
extern map<atom_t, size_t> atom_mu_loc;
extern vector<size_t> atom_mu_part_range;
extern int N_all_mu;

extern vector<atpair_t> tot_atpair;
extern vector<atpair_t> tot_atpair_ordered;
extern vector<atpair_t> local_atpair;

//! type alias of atom-pair mapping to real matrix indexed by unit-cell vector
typedef atom_mapping< map<Vector3_Order<int>, std::shared_ptr<matrix>> >::pair_t_old atpair_R_mat_t;
//! type alias of atom-pair mapping to complex matrix indexed by unit-cell vector
typedef atom_mapping< map<Vector3_Order<int>, std::shared_ptr<ComplexMatrix>> >::pair_t_old atpair_R_cplx_mat_t;
//! type alias of atom-pair mapping to complex matrix indexed by reciprocal vector
typedef atom_mapping< map<Vector3_Order<double>, std::shared_ptr<ComplexMatrix>> >::pair_t_old atpair_k_cplx_mat_t;

//! Tri-coefficient of localized RI (LRI) in real space. ABF on the first atom of the atom pair
extern atpair_R_mat_t Cs;
//! Coulomb matrix in ABF, represented in the reciprocal space.
extern atpair_k_cplx_mat_t Vq;
//! Truncated Coulomb matrix in ABF, represented in the reciprocal space.
extern atpair_k_cplx_mat_t Vq_cut;

extern map<Vector3_Order<double>, ComplexMatrix> Vq_block_loc;
extern map<Vector3_Order<double>, ComplexMatrix> Vq_cut_block_loc;

void allreduce_2D_coulomb_to_atompair(map<Vector3_Order<double>, ComplexMatrix> &Vq_loc, atpair_k_cplx_mat_t &coulomb_mat, double threshold );
void allreduce_atp_coulomb( atpair_k_cplx_mat_t &coulomb_mat );

int atom_iw_loc2glo(const int &atom_index, const int &iw_lcoal);
int atom_mu_loc2glo(const int &atom_index, const int &mu_lcoal);
int atom_mu_glo2loc(const int &glo_index, int &mu_index); //in-out mu_index, return atom_index;

void init_N_all_mu();
void allreduce_atp_aux();
//! inverse Fouriter transform of atom-pair mapping to complex matrix
atpair_R_cplx_mat_t inverse_FT_atpair_cplx_mat(atpair_k_cplx_mat_t kmat, vector<Vector3_Order<int>> Rlist);

vector<int> get_part_range();

//! Reshape Cs matrix from (n1*n2,n3) to (n2,n1*n3)
matrix reshape_Cs(size_t n1, size_t n2, size_t n3, const shared_ptr<matrix> &Csmat);
//! Reshape Cs matrix from (n1,n2*n3) to (n2,n1*n3)
matrix reshape_mat(size_t n1, size_t n2, size_t n3, const matrix &Csmat);
//! Reshape Cs matrix from n1,n2*n3) to (n1*n2,n3)
matrix reshape_mat_21(const size_t n1, const size_t n2, const size_t n3, const matrix &Csmat);
//! Reshape Cs matrix from n1,n2*n3) to (n1*n2,n3)
matrix reshape_mat_21(const size_t n1, const size_t n2, const size_t n3, const matrix &Csmat);
#endif
