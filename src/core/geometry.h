#pragma once
#include <array>

#include "atoms.h"

typedef std::array<double, 3> coord_t;

//! Cartesian coordinates of atoms
extern std::map<atom_t, coord_t> coord;

//! Fractional coordinates of atoms
extern std::map<atom_t, coord_t> coord_frac;

//! Atom type of atoms
extern std::map<atom_t, int> atom_types;
