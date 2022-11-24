/*!
 @file atoms.h
 @brief Utilies to handle atomic model and related data
 @date 2022-05-05 minyez created
 */
#ifndef ATOMS_H
#define ATOMS_H

#include <set>
#include <vector>
#include <utility>
#include <map>
/* #include "vector3.h" */
/* #include "matrix3.h" */

using std::vector;
using std::map;
using std::pair;
using std::size_t;

//! type of atom indices
typedef size_t atom_t;

//! atom-pair type. NOTE: may turn into a class?
typedef pair<atom_t, atom_t> atpair_t;

//! mapping from atom entries to data
template <typename T>
struct atom_mapping
{
    typedef T non_t;
    //! mapping between single atom and data
    typedef map<atom_t, T> single_t;
    //! mapping between atom pair and data
    typedef map<atpair_t, T> pair_t;
    //! mapping between tri-atom group and data
    typedef map<pair<atom_t, atpair_t>, T> tri_t;
    //! mapping between atom pair and data. Nested map, old style
    typedef map<atom_t, map<atom_t, T>> pair_t_old;
    //! mapping between tri-atom group and data. Nested map, old style
    typedef map<atom_t, map<atom_t, map<atom_t, T>>> tri_t_old;
};

//! get atom indcies from a atom_mapping::single_t object (old style)
template <typename T>
vector<atom_t> get_atom(const map<atom_t, T> &adata) // work
/* vector<atom_t> get_atom(const typename atom_mapping<T>::single_t &adata) // not work */
{
    vector<atom_t> alist;
    for (const auto &a1: adata)
            alist.push_back(a1.first);
    return alist;
}

//! get atom pairs from a atom_mapping::pair_t object (old style)
template <typename T>
vector<atpair_t> get_atom_pair(const map<atom_t, map<atom_t, T>> &apdata) // work
/* vector<pair<atom_t, atom_t>> get_atom_pair(const typename atom_mapping<T>::pair_t &apdata) // not work */
{
    vector<atpair_t> apair;
    for (const auto &a1: apdata)
        for (const auto &a2: a1.second)
            apair.push_back({a1.first, a2.first});
    return apair;
}

//! get atom pairs from a atom_mapping::pair_t object
template <typename T>
vector<atpair_t> get_atom_pair(const map<atpair_t, T> &apdata) // work
/* vector<pair<atom_t, atom_t>> get_atom_pair(const typename atom_mapping<T>::pair_t &apdata) // not work */
// FIXME: check this: https://en.cppreference.com/w/cpp/language/template_argument_deduction
{
    vector<atpair_t> apair;
    for (const auto &ap: apdata)
        apair.push_back(ap.first);
    return apair;
}

template <typename T>
vector<atpair_t> generate_atom_pair_from_nat(const T &nat, bool ordered_pair = false)
{
    vector<atpair_t> apair;
    for(int i=0; i!=nat; i++)
        for(int j = ordered_pair? 0 : i; j!=nat; j++)
            apair.push_back({i,j});
    return apair;
}

template <typename TA1, typename TA2>
std::pair<std::set<TA1>, std::set<TA2>> convert_IJset_to_Iset_Jset(std::set<std::pair<TA1, TA2>> IJset)
{
    std::set<TA1> Iset, Jset;
    for (const auto& IJ: IJset)
    {
        Iset.insert(IJ.first);
        Jset.insert(IJ.second);
    }
    return {Iset, Jset};
}

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
