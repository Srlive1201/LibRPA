#include "ri.h"
#include "atoms.h"

int n_irk_points;
int natom;
int ncell;
map<Vector3_Order<double>, double> irk_weight;
map<int, int> atom_nw;
map<int, int> atom_mu;

map<atom_t, map<atom_t, map<Vector3_Order<int>, std::shared_ptr<matrix>>>> Cs;
map<atom_t, map<atom_t, map<Vector3_Order<double>, std::shared_ptr<ComplexMatrix>>>> Vq;

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

