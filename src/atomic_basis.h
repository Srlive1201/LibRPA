/*!
 * @file atomic_basis.h
 * @brief Utilities for handling atomic basis functions
 */
#pragma once
#include <map>
#include <utility>
#include <vector>

#include "atoms.h"

namespace LIBRPA {

/*! @enum
 *
 *  Control for phase and order convention for basis functions in the same shell
 */
enum class AngularConvention
{
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
    std::vector<std::size_t> nbs_;
    std::vector<std::size_t> part_range_;
    std::map<atom_t, std::vector<std::size_t>> global_indices_;

    void initialize();
public:
    //! Total number of atoms
    std::size_t n_atoms;
    //! Total number of basis functions
    std::size_t nb_total;

    // Constructors
    AtomicBasis(): initialized_(false), nbs_(), part_range_(), n_atoms(0), nb_total(0) {};
    AtomicBasis(const std::vector<std::size_t>& nbs);
    AtomicBasis(const std::vector<int>& atoms,
                const std::map<int, std::size_t>& map_atom_nb);
    AtomicBasis(const std::map<std::size_t, std::size_t>& iatom_nbs);

    //! Set number of basis functions for each atom
    void set(const std::vector<std::size_t>& nbs);
    void set(const std::map<std::size_t, std::size_t>& iatom_nbs);

    //! Get the global index of a certain basis function of an atom
    std::size_t get_global_index(const int& i_atom, const std::size_t& i_loc_b) const;

    //! Get the size of submatrix corresponding to atom pair i and j
    std::size_t get_pair_matrix_size(const int& i_atom, const int& j_atom) const
    {
        return nbs_[i_atom] * nbs_[j_atom];
    };

    //! Get the global indices of all basis functions on an atom
    std::vector<std::size_t> get_global_indices(const int& i_atom) const;

    //! Get the index of atom on which a certain basis function is located
    int get_i_atom(const std::size_t& i_glo_b) const;

    //! Get the local indices of a basis function from its global index
    void get_local_index(const std::size_t& i_glo_b, int& i_atom, int& i_loc_b) const;
    int get_local_index(const std::size_t& i_glo_b, const int& i_atom) const;
    std::pair<int, int> get_local_index(const std::size_t& i_glo_b) const;

    std::size_t get_atom_nb(const int& i_atom) const { return nbs_[i_atom]; }
    std::vector<std::size_t> get_atom_nbs() const { return nbs_; }
    const std::vector<std::size_t>& get_part_range() const { return part_range_; }
    bool initialized() const { return initialized_; }
};

/*!
 * @brief Get the size of matrix block corresponding to an atom pair.
 *
 * @param  [in]  atbasis_r  AtomicBasis object for row
 * @param  [in]  at_r       Index of atomic for row
 * @param  [in]  atbasis_c  AtomicBasis object for column
 * @param  [in]  at_c       Index of atomic for column
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
 * @param  [in]  atbasis    AtomicBasis object for both row and columns
 * @param  [in]  at_r       Index of atomic for row
 * @param  [in]  at_c       Index of atomic for column
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
 *
 * @retval       indices    2D indices of elements in the request atom-pair blocks
 */
std::vector<std::pair<size_t, size_t>> get_2d_mat_indices_atpair(const AtomicBasis &atbasis_r,
                                                                 const AtomicBasis &atbasis_c,
                                                                 const std::vector<atpair_t> &IJs,
                                                                 bool row_fast);

/*!
 * @brief Get the 1D indices for matrix elements in the atom-pair blocks
 *        corresponding to requested atom pairs.
 *
 * @param  [in]  atbasis_r  AtomicBasis object for row
 * @param  [in]  atbasis_c  AtomicBasis object for column
 * @param  [in]  IJs        Atomic pairs to request
 * @param  [in]  row_fast   Flag to set the row basis index faster
 * @param  [in]  row_major  Flag to compute 1D indices in row major (C-style)
 *
 * @retval       indices    1D indices of elements
 */
std::vector<size_t> get_1d_mat_indices_atpair(const AtomicBasis &atbasis_r,
                                              const AtomicBasis &atbasis_c,
                                              const std::vector<atpair_t> &IJs,
                                              bool row_fast,
                                              bool row_major);

extern AtomicBasis atomic_basis_wfc;
extern AtomicBasis atomic_basis_abf;

} // namespace LIBRPA
