#include "ri.h"

int n_irk_points;
int natom;
int ncell;
map<Vector3_Order<double>, double> irk_weight;
map<atom_t, size_t> atom_nw;
map<atom_t, size_t> atom_mu;

apair_R_mat_t Cs;
apair_k_cplx_mat_t Vq;

int atom_iw_loc2glo(const int &atom_index, const int &iw_lcoal)
{
    int nb = 0;
    for (int ia = 0; ia != atom_index; ia++)
        nb += atom_nw[ia];
    return iw_lcoal + nb;
}

int atom_mu_loc2glo(const int &atom_index, const int &mu_lcoal)
{
    int nb = 0;
    for (int ia = 0; ia != atom_index; ia++)
        nb += atom_mu[ia];
    return mu_lcoal + nb;
}

