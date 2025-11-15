#include "atomic_basis.h"

#include <algorithm>
#include <cassert>
#include <numeric>
#include <stdexcept>

#include "base_utility.h"

namespace LIBRPA {

void AtomicBasis::initialize()
{
    n_atoms = nbs_.size();

    if (n_atoms < 0)
        throw std::invalid_argument("empty basis is parsed");
    part_range_.resize(n_atoms + 1);
    part_range_[0] = 0;
    for (size_t i = 0; i < n_atoms; i++)
        part_range_[i+1] = part_range_[i] + nbs_[i];
    nb_total = part_range_[n_atoms];

    // compute local indices
    glo2iat_.resize(nb_total);
    glo2loc_.resize(nb_total);
    for (int I = 0; I < as_int(n_atoms); I++)
    {
        std::fill(glo2iat_.begin()+part_range_[I], glo2iat_.begin()+part_range_[I+1], I);
        std::iota(glo2loc_.begin()+part_range_[I], glo2loc_.begin()+part_range_[I+1], 0);
    }

    initialized_ = true;
}

AtomicBasis::AtomicBasis(const std::vector<int>& atoms,
                         const std::map<int, std::size_t>& map_atom_nb)
    : initialized_(false), nbs_(), part_range_(), glo2iat_(), glo2loc_()
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
    : initialized_(false), nbs_(nbs), part_range_(), n_atoms(0), nb_total(0)
{
    initialize();
}

AtomicBasis::AtomicBasis(const std::map<size_t, std::size_t>& iatom_nbs)
    : initialized_(false), nbs_(), part_range_(), n_atoms(0), nb_total(0)
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

std::vector<std::size_t> AtomicBasis::get_global_indices(const int& i_atom) const
{
    std::vector<std::size_t> indices(this->nbs_[i_atom]);
    std::iota(indices.begin(), indices.end(), get_part_range()[i_atom]);
    return indices;
}

std::vector<std::pair<size_t, size_t>>
get_2d_mat_indices_atpair(const AtomicBasis &atbasis_r,
                          const AtomicBasis &atbasis_c,
                          const std::vector<atpair_t> &IJs,
                          const bool row_fast,
                          const bool sort_fast)
{
    std::vector<std::pair<size_t, size_t>> indices;
    size_t nelems = 0;

    // pre-compute the size of all indices to reserve enough space
    for (const auto &IJ: IJs)
    {
        nelems += get_pair_matrix_size(atbasis_r, IJ.first, atbasis_c, IJ.second);
    }
    indices.reserve(nelems);

    const auto append_indices = [&](const atpair_t& IJ) {
        const auto row_ids = atbasis_r.get_global_indices(IJ.first);
        const auto col_ids = atbasis_c.get_global_indices(IJ.second);

        if (row_fast)
        {
            for (const auto &c : col_ids)
                for (const auto &r : row_ids) indices.push_back({r, c});
        }
        else
        {
            for (const auto &r : row_ids)
                for (const auto &c : col_ids) indices.push_back({r, c});
        }
    };

    if (sort_fast)
    {
        // pre-sort the atom pairs to make index sorting easier
        std::vector<atpair_t> IJs_sorted = IJs;
        std::sort(IJs_sorted.begin(), IJs_sorted.end(), FastLess<atpair_t>{row_fast});

        for (const auto &IJ: IJs_sorted) append_indices(IJ);
        std::sort(indices.begin(), indices.end(), FastLess<gloid_ap_t>{row_fast});
    }
    else
    {
        for (const auto &IJ: IJs) append_indices(IJ);
    }

    return indices;
}

std::vector<size_t> get_1d_mat_indices_atpair(const AtomicBasis &atbasis_r,
                                              const AtomicBasis &atbasis_c,
                                              const std::vector<atpair_t> &IJs,
                                              const bool row_fast,
                                              const bool row_major,
                                              const bool sort_fast)
{
    const auto indices_2d = get_2d_mat_indices_atpair(atbasis_r, atbasis_c, IJs, row_fast, sort_fast);
    const auto indices_1d = flatten_2d_indices(indices_2d, atbasis_r.nb_total, atbasis_c.nb_total, row_major);
    return indices_1d;
}

AtomicBasis atomic_basis_wfc;
AtomicBasis atomic_basis_abf;

} // namespace LIBRPA
