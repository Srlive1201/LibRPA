/*!
 * @file atomic_basis.h
 * @brief Utilities for handling atomic basis functions
 */
#pragma once
#include <cassert>
#include <map>
#include <utility>
#include <vector>

#include "atom.h"

namespace librpa_int {

typedef int locid_t;
typedef std::size_t gloid_t;
typedef std::pair<locid_t, locid_t> locid_ap_t;
typedef std::pair<gloid_t, gloid_t> gloid_ap_t;

/*! @enum
 *
 *  Control for phase and order convention for basis functions in the same shell
 */
enum class AngularConvention
{
    UNSET = -1,
    AIMS = 0,
    ABACUS = 1,
    OPENMX = 2,
    PYSCF = 3,
};

/*! @class
 * @brief Object to handle atomic basis
 */
class AtomicBasis
{
private:
    bool initialized_;
public:
    std::vector<std::size_t> nbs_;
private:
    std::vector<std::size_t> part_range_;
    std::vector<int> glo2iat_;
    std::vector<std::size_t> glo2loc_;
    // std::vector<int> irad_;
    // std::vector<int> l_;
    // std::vector<int> m_;

    void initialize();
public:
    //! Total number of atoms
    std::size_t n_atoms;
    //! Total number of basis functions
    std::size_t nb_total;

    // Constructors
    AtomicBasis(): initialized_(false), nbs_(), part_range_(), glo2iat_(), glo2loc_(), n_atoms(0), nb_total(0) {};
    AtomicBasis(const std::vector<std::size_t>& nbs);
    AtomicBasis(const std::vector<int>& atom_species,
                const std::map<int, std::size_t>& map_species_nb);
    AtomicBasis(const std::map<std::size_t, std::size_t>& iatom_nbs);

    //! Set number of basis functions for each atom
    void set(const std::vector<std::size_t>& nbs);
    void set(const std::vector<int>& atom_species,
             const std::map<int, std::size_t>& map_species_nb);
    void set(const std::map<std::size_t, std::size_t>& iatom_nbs);

    //! Get the global index of a certain basis function of an atom
    inline std::size_t get_global_index(const int& i_atom, const std::size_t& i_loc_b) const noexcept
    {
        return part_range_[i_atom] + i_loc_b;
    }

    //! Get the size of submatrix corresponding to atom pair i and j
    std::size_t get_pair_matrix_size(const int& i_atom, const int& j_atom) const noexcept
    {
        return nbs_[i_atom] * nbs_[j_atom];
    };

    //! Get the global indices of all basis functions on an atom
    std::vector<std::size_t> get_global_indices(const int& i_atom) const;

    //! Get the index of atom on which a certain basis function is located
    inline int get_i_atom(const std::size_t& i_glo_b) const noexcept
    {
        assert (i_glo_b < nb_total && i_glo_b >= 0);
        return glo2iat_[i_glo_b];
    }

    //! Get the local indices of a basis function from its global index
    inline void get_local_index(const std::size_t& i_glo_b, int& i_atom, int& i_loc_b) const
    {
        i_atom = get_i_atom(i_glo_b);
        i_loc_b = as_int(glo2loc_[i_glo_b]);
    }
    inline void get_local_index(const std::size_t& i_glo_b, size_t& i_atom, size_t& i_loc_b) const
    {
        i_atom = as_size(get_i_atom(i_glo_b));
        i_loc_b = glo2loc_[i_glo_b];
    }
    inline int get_local_index(const std::size_t& i_glo_b, const int& i_atom) const noexcept
    {
        // return i_glo_b - part_range_[i_atom];
        return as_int(glo2loc_[i_glo_b]);
    }

    inline std::pair<int, int> get_local_index(const std::size_t& i_glo_b) const noexcept
    {
        // int i_atom, i_loc_b;
        // this->get_local_index(i_glo_b, i_atom, i_loc_b);
        return {glo2iat_[i_glo_b], as_int(glo2loc_[i_glo_b])};
    }

    inline std::size_t get_atom_nb(int i_atom) const { return nbs_.at(as_size(i_atom)); }
    inline std::size_t operator[](int i_atom) const { return nbs_.at(as_size(i_atom)); }
    inline std::vector<std::size_t> get_atom_nbs() const noexcept { return nbs_; }
    inline const std::vector<std::size_t>& get_part_range() const noexcept { return part_range_; }
    inline bool initialized() const noexcept { return initialized_; }
};

/*!
 * @brief Get the size of matrix block corresponding to an atom pair.
 *
 * @param  [in]  atbasis_r  AtomicBasis object for row
 * @param  [in]  at_r       Index of atom for row
 * @param  [in]  atbasis_c  AtomicBasis object for column
 * @param  [in]  at_c       Index of atom for column
 *
 * @retval       size       Number of elements in the atom-pair matrix block.
 */
inline std::size_t get_pair_matrix_size(const AtomicBasis& atbasis_r, const int& at_r,
                                        const AtomicBasis& atbasis_c, const int& at_c)
{
    return atbasis_r.get_atom_nb(at_r) * atbasis_c.get_atom_nb(at_c);
}

/*!
 * @brief Get the size of matrix block corresponding to an atom pair.
 *
 * @param  [in]  atbasis    AtomicBasis object for both row and column
 * @param  [in]  at_r       Index of atom for row
 * @param  [in]  at_c       Index of atom for column
 *
 * @retval       size       Number of elements in the atom-pair matrix block.
 */
inline std::size_t get_pair_matrix_size(const AtomicBasis& atbasis, const int& at_r, const int& at_c)
{
    return get_pair_matrix_size(atbasis, at_r, atbasis, at_c);
}

/*!
 * @brief Get the 2D indices for matrix elements in the atom-pair blocks
 *        corresponding to requested atom pairs.
 *
 * @param  [in]  atbasis_r  AtomicBasis object for row
 * @param  [in]  atbasis_c  AtomicBasis object for column
 * @param  [in]  IJs        Atomic pairs to request
 * @param  [in]  row_fast   Flag to set the row basis index faster
 * @param  [in]  sort_fast  Flag to sort with the faster dimension.
 *                          The indices are continuous within one atom-pair block.
 *                          When sort_fast is true, the indices will be sorted and
 *                          hence the continuity may break.
 *
 * @retval       indices    2D indices of elements in the request atom-pair blocks
 */
std::vector<std::pair<size_t, size_t>> get_2d_mat_indices_atpair(const AtomicBasis &atbasis_r,
                                                                 const AtomicBasis &atbasis_c,
                                                                 const std::vector<atpair_t> &IJs,
                                                                 const bool row_fast,
                                                                 const bool sort_fast = false);

/*!
 * @brief Get the 1D indices for matrix elements in the atom-pair blocks
 *        corresponding to requested atom pairs.
 *
 * @param  [in]  atbasis_r  AtomicBasis object for row
 * @param  [in]  atbasis_c  AtomicBasis object for column
 * @param  [in]  IJs        Atomic pairs to request
 * @param  [in]  row_fast   Flag to set the row basis index faster
 * @param  [in]  row_major  Flag to compute 1D indices in row major (C-style)
 * @param  [in]  sort_fast  Flag to sort with the faster dimension.
 *                          The indices are continuous within one atom-pair block.
 *                          When sort_fast is true, the indices will be sorted and
 *                          hence the continuity may break.
 *
 * @retval       indices    1D indices of elements
 */
std::vector<size_t> get_1d_mat_indices_atpair(const AtomicBasis &atbasis_r,
                                              const AtomicBasis &atbasis_c,
                                              const std::vector<atpair_t> &IJs,
                                              const bool row_fast,
                                              const bool row_major,
                                              const bool sort_fast = false);

// extern AtomicBasis atomic_basis_wfc;
// extern AtomicBasis atomic_basis_abf;

} // namespace librpa_int
