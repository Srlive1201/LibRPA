#include <cassert>
#include <stdexcept>
#include <string>
#include "atomic_basis.h"

namespace LIBRPA {

void AtomicBasis::initialize()
{
    n_atoms = nbs_.size();

    if (n_atoms < 0)
        throw std::invalid_argument("empty basis is parsed");
    part_range.resize(n_atoms);
    part_range[0] = 0;
    for (int i = 0; i != n_atoms - 1; i++)
        part_range[i+1] = part_range[i] + nbs_[i];
    nb_total = part_range[n_atoms-1] + nbs_[n_atoms-1];
}

AtomicBasis::AtomicBasis(const std::vector<int>& atoms,
                         const std::map<int, std::size_t>& map_atom_nb)
    : nbs_()
{
    for (const auto &atom: atoms)
    {
        if (map_atom_nb.count(atom))
            this->nbs_.push_back(map_atom_nb.at(atom));
    }
    if (atoms.size() != this->nbs_.size())
    {
        throw std::invalid_argument("missing nbs for some atom");
    }
    initialize();
}

AtomicBasis::AtomicBasis(const std::vector<std::size_t>& nbs)
    : nbs_(nbs)
{
    initialize();
}

void AtomicBasis::set(const std::vector<std::size_t>& nbs)
{
    nbs_.clear();
    nbs_ = nbs;
    initialize();
}

std::size_t AtomicBasis::get_global_index(const int& i_atom, const std::size_t& i_loc_b) const
{
    return part_range[i_atom] + i_loc_b;
}

int AtomicBasis::get_i_atom(const std::size_t& i_glo_b) const
{
    if (i_glo_b >= nb_total)
        throw std::invalid_argument("Global index out-of-bound: " + std::to_string(i_glo_b));
    int iat;
    for (iat = 0; iat < n_atoms; iat++)
        if (i_glo_b < (part_range[iat] + nbs_[iat]))
            break;
    return iat;
}

void AtomicBasis::get_local_index(const std::size_t& i_glo_b, int& i_atom, int& i_loc_b) const
{
    i_atom = get_i_atom(i_glo_b);
    i_loc_b = i_glo_b - part_range[i_atom];
}

std::pair<int, int> AtomicBasis::get_local_index(const std::size_t& i_glo_b) const
{
    int i_atom, i_loc_b;
    i_atom = get_i_atom(i_glo_b);
    i_loc_b = i_glo_b - part_range[i_atom];
    return {i_atom, i_loc_b};
}

AtomicBasis atomic_basis_wfc;
AtomicBasis atomic_basis_abf;

} // namespace LIBRPA
