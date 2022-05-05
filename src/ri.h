#ifndef RI_H
#define RI_H

#include <map>
#include <memory>
#include <utility>
#include "vector3_order.h"
#include "atoms.h"
#include "complexmatrix.h"
using std::map;
// using std::pair; // TODO: use pair as first map

extern int n_irk_points;
extern int natom;
extern int ncell;
extern map<Vector3_Order<double>, double> irk_weight;
extern map<atom_t, size_t> atom_nw;
extern map<atom_t, size_t> atom_mu;

extern map<atom_t, map<atom_t, map<Vector3_Order<int>, std::shared_ptr<matrix>>>> Cs;
extern map<atom_t, map<atom_t, map<Vector3_Order<double>, std::shared_ptr<ComplexMatrix>>>> Vq;

int atom_iw_loc2glo(const int &atom_index, const int &iw_lcoal);
int atom_mu_loc2glo(const int &atom_index, const int &mu_lcoal);

#endif // !RI_H
