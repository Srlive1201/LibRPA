/*!
 @file atoms.h
 @brief Utilies to handle atomic model and related data
 @date 2022-05-05 minyez created
 */
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

//! mapping from atom entries to data
template <typename T>
struct atom_mapping
{
    //! mapping between single atom and data
    typedef map<atom_t, T> single_t;
    //! mapping between atom pair and data
    typedef map<pair<atom_t, atom_t>, T> pair_t;
    //! mapping between tri-atom group and data
    typedef map<pair<atom_t, pair<atom_t, atom_t>>, T> tri_t;
};

//! Object to handle the atomic structure
/*!
  It should also handle the non-PBC case.
  @note For now, atom positions are not necessary, as one decides the sparsity from matrix element,
  instead of distances between atoms
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
