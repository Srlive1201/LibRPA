#include "geometry.h"
#include <stdexcept>

namespace librpa_int {

std::map<atom_t, coord_t> coord;
std::map<atom_t, coord_t> coord_frac;

std::map<atom_t, int> atom_types;

void Atoms::set(const std::vector<int> &types_in,
                const std::vector<coord_t> &coords_in)
{
    auto n_type_orig = this->types.size();
    auto n_coord_orig = this->coords.size();
    auto n_type_in = types_in.size();
    auto n_coord_in = coords_in.size();

    // Do nothing
    if (n_coord_in == 0 && n_type_in == 0) return;

    auto set_types = [&]()
    {
        this->types.clear();
        for (size_t iat = 0; iat < n_type_in; iat++)
        {
            this->types[atom_t(iat)] = types_in[iat];
        }
    };
    auto set_coords = [&]()
    {
        this->coords.clear();
        this->coords_frac.clear();
        for (size_t iat = 0; iat < n_coord_in; iat++)
        {
            this->coords[atom_t(iat)] = coords_in[iat];
        }
    };

    auto set_both = [&]()
    {
        set_coords();
        set_types();
    };
    auto set_coords_with_known_types = [&]()
    {
        if (n_coord_in != n_type_orig)
            throw std::runtime_error("Number of input coords inconsistent with parsed types");
        set_coords();
    };
    auto set_types_with_known_coords = [&]()
    {
        if (n_coord_orig != n_type_in)
            throw std::runtime_error("Number of input types inconsistent with parsed coords");
        set_types();
    };

    // When both inputs are non-empty: check consistency and set both
    if (n_type_in > 0 && n_coord_in > 0)
    {
        if (n_type_in != n_coord_in) throw std::runtime_error("Imbalanced coords and types inputs");
        set_both(); return;
    }

    // Structure is not set yet: set anyway
    if (n_type_orig == 0 && n_coord_orig == 0)
    {
        set_both(); return;
    }

    // Coordinates set, but types were not: either reset coords, or set types
    if (n_type_orig == 0)
    {
        n_type_in > 0 ? set_types_with_known_coords() : set_coords();
        return;
    }

    // Types set, but coordinates were not: Either reset types, or set coords
    if (n_coord_orig == 0)
    {
        n_coord_in > 0 ? set_coords_with_known_types() : set_types();
        return;
    }

    // Both are already set: reset coords/types if parsed
    n_coord_in > 0 ? set_coords_with_known_types() : set_types_with_known_coords();
}

void Atoms::set(const std::vector<int> &types_in,
                const std::vector<coord_t> &coords_in,
                const Matrix3 &lattice)
{
    this->set(types_in, coords_in);
    if (!is_set())
    {
        throw std::runtime_error("Coordinates not set, cannot compute fractional coords");
    }

    const auto G = lattice.Inverse().Transpose();

    this->coords_frac.clear();
    auto n_atoms = this->coords.size();
    for (size_t i = 0; i < n_atoms; i++)
    {
        this->coords_frac[atom_t(i)] = {
            (this->coords[i][0] * G.e11 + this->coords[i][1] * G.e12 + this->coords[i][2] * G.e13),
            (this->coords[i][0] * G.e21 + this->coords[i][1] * G.e22 + this->coords[i][2] * G.e23),
            (this->coords[i][0] * G.e31 + this->coords[i][1] * G.e32 + this->coords[i][2] * G.e33),
        };
    }
}


}
