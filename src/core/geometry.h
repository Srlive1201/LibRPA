#pragma once
#include <array>

#include "../math/matrix3.h"
#include "atom.h"

namespace librpa_int {

typedef std::array<double, 3> coord_t;

//! Cartesian coordinates of atoms
extern std::map<atom_t, coord_t> coord;

//! Fractional coordinates of atoms
extern std::map<atom_t, coord_t> coord_frac;

//! Atom type of atoms
extern std::map<atom_t, int> atom_types;

//! Object to handle the atomic structure
class Atoms
{
private:
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

    void set(const std::vector<int> &types_in,
             const std::vector<coord_t> &coords_in);

    void set(const std::vector<int> &types_in,
             const std::vector<coord_t> &coords_in,
             const Matrix3 &lattice);

    bool is_set() { return coords.size() > 0; }
    bool is_type_set() { return types.size() > 0; }
    bool is_frac_set() { return coords_frac.size() > 0; }
};

}
