#include "atomic_basis.h"

#include <algorithm>
#include <cassert>
#include <functional>
#include <stdexcept>
#include <string>

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
    : nbs_(), part_range(), n_atoms(0), nb_total(0)
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
    : nbs_(nbs), part_range(), n_atoms(0), nb_total(0)
{
    initialize();
}

AtomicBasis::AtomicBasis(const std::map<size_t, std::size_t>& iatom_nbs)
    : nbs_(), part_range(), n_atoms(0), nb_total(0)
{
    // sort atom index first
    // std::function<bool(const std::pair<std::size_t, std::size_t>&,
    //                    const std::pair<std::size_t, std::size_t>&)>
    set(iatom_nbs);
}

void AtomicBasis::set(const std::vector<std::size_t>& nbs)
{
    nbs_.clear();
    nbs_ = nbs;
    initialize();
}

void AtomicBasis::set(const std::map<std::size_t, std::size_t>& iatom_nbs)
{
    auto sortfn = [](const std::pair<std::size_t, std::size_t> &p1,
                     const std::pair<std::size_t, std::size_t>& p2)
                    { return p1.first < p2.first; };
    std::vector<std::pair<std::size_t, std::size_t>> v_ianb;
    for (const auto &ia_nb: iatom_nbs)
        v_ianb.push_back(ia_nb);
    std::sort(v_ianb.begin(), v_ianb.end(), sortfn);
    std::vector<std::size_t> nbs;
    nbs_.clear();
    for (const auto &nb: v_ianb)
        nbs_.push_back(nb.second);
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

int AtomicBasis::get_local_index(const std::size_t& i_glo_b, const int& i_atom) const
{
    return i_glo_b - part_range[i_atom];
}

std::pair<int, int> AtomicBasis::get_local_index(const std::size_t& i_glo_b) const
{
    int i_atom, i_loc_b;
    this->get_local_index(i_glo_b, i_atom, i_loc_b);
    return {i_atom, i_loc_b};
}

AtomicBasis atomic_basis_wfc;
AtomicBasis atomic_basis_abf;

} // namespace LIBRPA
