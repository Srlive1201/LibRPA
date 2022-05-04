#ifndef ATOMS_H
#define ATOMS_H

#include <vector>
#include <utility>
#include <map>
/* #include "vector3.h" */
/* #include "matrix3.h" */

/* using std::vector; */
using std::map;
using std::pair;
using std::size_t;

//! type of atom indices
typedef size_t atom_t;

//! a general mapping from atom pairs to any type
template <typename T>
struct atompair_map
{
    typedef map<pair<atom_t, atom_t>, T> type;
};

//! Object to handle the atomic structure
/*!
  It should also handle the non-PBC case.
  For now atom positions are not necessary, as one decides the sparsity from matrix element,
  instead of distances.
 */
class Atoms
{
    private:
        //! number of inequivalent atoms
        /* size_t n_ineq; */
        //! total number of atoms
        /* Matrix3 latt; */
        /* Matrix3 G; */
        size_t n_atoms;
        //! number of atoms of each inequivalent type;
        /* vector<size_t> n_atoms_per_ineq; */
        //! cooridnates of all atoms, Cartesian
        /* vector<Vector3<double>> coords; */
        //! cooridnates of the representative inequivalent atoms, Cartesian
        /* vector<Vector3<double>> coords_ineq; */
        //! index mapping from total atom to inquivalent type
        /* vector<map<atom_t, atom_t>> imap_at2ineq; */
    public:
        /* Matrix3 & get_lattice_vectors() { return latt; } */
        /* const Matrix3 & get_lattice_vectors() const { return latt; } */
        /* const size_t & get_n_atoms() const { return n_atoms; } */
};

#endif // !ATOMS_H
