#pragma once
#include <array>

#include "../math/matrix3.h"
#include "atom.h"

namespace librpa_int {

typedef std::array<double, 3> coord_t;

// //! Cartesian coordinates of atoms
// extern std::map<atom_t, coord_t> coord;

// //! Fractional coordinates of atoms
// extern std::map<atom_t, coord_t> coord_frac;

// //! Atom type of atoms
// extern std::map<atom_t, int> atom_types;

//! Object to handle the atomic structure
class Atoms
{
private:
    static Matrix3 pseudo_recplatt_in_2pi_;

private:
    bool is_pseudo_lattice_ = true;
    // //! number of atoms of each inequivalent type;
    // /* vector<size_t> n_atoms_per_ineq; */
    // //! cooridnates of all atoms, Cartesian
    // //! cooridnates of the representative inequivalent atoms, Cartesian
    // /* vector<Vector3<double>> coords_ineq; */
    // //! index mapping from total atom to inquivalent type
    // /* vector<map<atom_t, atom_t>> imap_at2ineq; */

public:
    std::map<atom_t, int> types;
    std::map<atom_t, coord_t> coords;
    std::map<atom_t, coord_t> coords_frac;

    Atoms() = default;
    Atoms(const std::vector<int> &types_in, const std::vector<coord_t> &coords_in) { set(types_in, coords_in); }
    ~Atoms() = default;

    size_t size() const { return std::max(types.size(), coords.size()); }

    void set(const std::vector<int> &types_in,
             const std::vector<coord_t> &coords_in);

    void set(const std::vector<int> &types_in,
             const std::vector<coord_t> &coords_in,
             const Matrix3 &lattice);

    bool is_set() const { return coords.size() > 0; }
    bool is_type_set() const { return types.size() > 0; }
    bool is_frac_set() const { return !is_pseudo_lattice_; }
};

}
