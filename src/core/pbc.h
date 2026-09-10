/*!
 @file pbc.h
 @brief Utilities to deal with periodic boundary conditions
 */
#pragma once
#include <array>
#include <cstddef>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

#include "../math/matrix3.h"
#include "../math/vector3_order.h"

namespace librpa_int {

//! Relation between one stored full-grid k-point and the loaded SCF k-list.
enum class KFullToKRelation
{
    //! No loaded SCF k-point represents this full-grid point yet.
    NONE,
    //! The loaded SCF k-point is the same k-point.
    DIRECT,
    //! The loaded SCF k-point is the time-reversal partner, -k.
    TIME_REVERSAL,
};

class PeriodicBoundaryData
{
private:
    bool lattice_reset_ = false;
    bool kgrid_no_symmetry_ = true;

    void reset_k_info();

public:
    //! Lattice vectors as a 3D-matrix, each row as a lattice vector. Unit: Bohr
    Matrix3 latvec;
    //! Reciprocal lattice vectors as a 3D-matrix, each row as a reciprocal vector. Unit: 2pi/Bohr
    Matrix3 G;
    //! Lattice vectors as a nested array for LibRI call
    std::array<std::array<double, 3>, 3> latvec_array;

    //! Period of each lattice vector in the BvK cell
    Vector3_Order<int> period;
    //! Period as an array for LibRI call
    std::array<int, 3> period_array;

    std::vector<Vector3_Order<int>> Rlist;
    //! Loaded SCF k-points in Cartesian coordinates
    std::vector<Vector3_Order<double>> klist;
    //! Loaded SCF k-points in fractional coordinates
    std::vector<Vector3_Order<double>> kfrac_list;
    //! Scalar weight of each loaded SCF k-point, if available.
    std::vector<double> weight_k;
    //! Full-BZ k-points when explicitly expanded; otherwise mirrors klist.
    std::vector<Vector3_Order<double>> klist_full;
    //! Full-BZ k-points in fractional coordinates when explicitly expanded; otherwise mirrors kfrac_list.
    std::vector<Vector3_Order<double>> kfrac_list_full;
    //! Mapping from loaded SCF k-points to stored full-grid k-points.
    std::vector<int> k_to_kfull;
    //! Mapping from stored full-grid k-points to loaded SCF representative k-points; -1 means missing.
    std::vector<int> kfull_to_k;
    //! How each kfull_to_k entry was resolved.
    std::vector<KFullToKRelation> kfull_to_k_relation;
    //! Whether any missing full-grid k-point was found through time reversal.
    bool kgrid_uses_time_reversal = false;

    //! K-points where Coulomb matrices are parsed, in Cartesian coordinates
    std::vector<Vector3_Order<double>> klist_coul;
    //! Scalar weight of each Coulomb k-point
    std::vector<double> weight_q;
    //! Same as weight_q, but indexed directly by the Coulomb k-point
    std::map<Vector3_Order<double>, double> map_q_weight;
    std::vector<int> irk_point_id_mapping;
    std::vector<int> isymops;
    //! Mapping of Coulomb k-points to loaded SCF k-points using that representative
    std::map<Vector3_Order<double>, std::vector<Vector3_Order<double>>> map_irk_ks;

    PeriodicBoundaryData();

    //! Setting lattice vectors. latt_mat should be in Bohr unit.
    void set_latvec(const std::vector<double> &lat_mat);
    //! Set lattice and reciprocal lattice vectors. The parsed rec_mat should be in Bohr^-1 unit. No consistency check.
    void set_latvec_and_G(const std::vector<double> &latt_mat,
                          const std::vector<double> &recp_mat);
    //! Set BvK periodicity and rebuild the corresponding R-grid.
    void set_period(int nk1, int nk2, int nk3);
    //! Set the loaded SCF k-point list. kvecs should be in Bohr^-1 unit.
    void set_kgrids_kvec(int nk1, int nk2, int nk3,
                         const std::vector<double> &kvecs,
                         const std::vector<double> &kweights = {});
    //! Set an irreducible loaded k-point list plus full-BZ star members in internal units.
    void set_irreducible_kgrids_kvec(
        int nk1, int nk2, int nk3,
        const std::vector<double> &kvecs_ibz,
        const std::vector<std::vector<Vector3_Order<double>>> &full_kstars);
    //! Set the mapping from loaded SCF k-points to Coulomb q-points.
    void set_kq_mapping(const std::vector<int> &map_q_ks,
                         const std::vector<int> &isymops_in = {});

    // Getting
    int get_R_index(const Vector3_Order<int> &R) const;
    int get_k_index_full(const Vector3_Order<double> &k) const;
    int get_k_index_ibz(const Vector3_Order<double> &k) const;
    int get_n_cells_bvk() const { return period.x * period.y * period.z; }
    bool is_latt_set() const { return lattice_reset_; }
};

// TODO: make it into a template
std::vector<Vector3_Order<int>> construct_R_grid(const Vector3_Order<int> &period);

std::vector<Vector3_Order<double>> build_uniform_kmesh_frac(const Vector3_Order<int> &period);

//! Get the index of R in an Rlist. If R is not found in the list, return a negative number
int get_R_index(const std::vector<Vector3_Order<int>> &Rlist, const Vector3_Order<int> &R);

bool is_gamma_point(const Vector3_Order<double> &kpt, double thres = 1.0e-5);
bool is_gamma_point(const Vector3_Order<int> &kpt);

Vector3_Order<int> find_nearest_bvk_cell(const Vector3<double> &coord_frac_I,
                                         const Vector3<double> &coord_frac_J,
                                         const Vector3_Order<int> &bvk_direct,
                                         const Vector3_Order<int> &period, const Matrix3 &latvec);

std::vector<Vector3_Order<int>> find_nearest_bvk_cells(const Vector3<double> &coord_frac_I,
                                                       const Vector3<double> &coord_frac_J,
                                                       const Vector3_Order<int> &bvk_direct,
                                                       const Vector3_Order<int> &period,
                                                       const Matrix3 &latvec);

// Remap parsed lattice vectors to their BvK counterparts for each atom pair.
// remap_convention = 0: choose one nearest BvK cell.
// remap_convention = 1: keep all nearest BvK cells sharing the same minimal distance.
template <typename AtomT, typename AtomPairT = std::pair<AtomT, AtomT>>
class AtomPairBvKRemap
{
public:
    using atom_type = AtomT;
    using atom_pair_type = AtomPairT;
    using R_type = Vector3_Order<int>;
    using R_remap_type = std::vector<R_type>;
    using remap_type = std::map<atom_pair_type, std::map<R_type, R_remap_type>>;
    using const_iterator = typename remap_type::const_iterator;

    AtomPairBvKRemap() = default;

    template <typename CoordT>
    AtomPairBvKRemap(const std::map<atom_type, CoordT> &coord_fracs,
                     const std::vector<R_type> &Rs,
                     const R_type &period, const Matrix3 &latvec,
                     int remap_convention = 0)
    {
        this->build(coord_fracs, Rs, period, latvec, remap_convention);
    }

    template <typename CoordT>
    void build(const std::map<atom_type, CoordT> &coord_fracs,
               const std::vector<R_type> &Rs,
               const R_type &period, const Matrix3 &latvec,
               int remap_convention = 0)
    {
        remap_.clear();
        if (remap_convention != 0 && remap_convention != 1)
            throw std::runtime_error("Invalid BvK remap convention");

        for (const auto &[I, coord_frac_I]: coord_fracs)
        {
            for (const auto &[J, coord_frac_J]: coord_fracs)
            {
                for (const auto &R: Rs)
                {
                    R_remap_type R_bvks;
                    if (remap_convention == 0)
                    {
                        R_bvks.push_back(find_nearest_bvk_cell(coord_frac_I, coord_frac_J, R, period, latvec));
                    }
                    else
                    {
                        R_bvks = find_nearest_bvk_cells(coord_frac_I, coord_frac_J, R, period, latvec);
                    }

                    if (R_bvks.size() == 1 && R_bvks.front() == R) continue;
                    remap_[atom_pair_type{I, J}][R] = std::move(R_bvks);
                }
            }
        }
    }

    const remap_type &data() const { return remap_; }
    const R_remap_type *find_R_bvk(const atom_pair_type &atom_pair, const R_type &R) const
    {
        const auto it_atom_pair = remap_.find(atom_pair);
        if (it_atom_pair == remap_.cend()) return nullptr;

        const auto it_R = it_atom_pair->second.find(R);
        if (it_R == it_atom_pair->second.cend()) return nullptr;

        return &it_R->second;
    }
    const std::map<R_type, R_remap_type> &at(const atom_pair_type &atom_pair) const { return remap_.at(atom_pair); }
    std::size_t size() const { return remap_.size(); }
    bool empty() const { return remap_.empty(); }
    bool is_empty() const { return remap_.empty(); }
    void clear() { remap_.clear(); }
    const_iterator begin() const { return remap_.begin(); }
    const_iterator end() const { return remap_.end(); }

private:
    remap_type remap_;
};

// extern int kv_nmp[3];
// //! lattice vectors as a 3D-matrix, each row as a lattice vector. Unit: Bohr
// extern Matrix3 latvec;
// //! same as latvec, but a nested array for LibRI call
// extern std::array<std::array<double, 3>, 3> lat_array;
// //! reciprocal lattice vectors as a 3D-matrix, each row as a reciprocal vector. Unit: 2pi/Bohr
// extern Matrix3 G;
// extern std::vector<Vector3_Order<double>> klist;
// extern std::vector<Vector3_Order<double>> klist_coul;
// extern std::vector<Vector3_Order<double>> kfrac_list;
// extern std::vector<int> irk_point_id_mapping;
// extern std::map<Vector3_Order<double>, std::vector<Vector3_Order<double>>> map_irk_ks;
// extern Vector3<double> *kvec_c;

int get_k_index_full(const Vector3_Order<double> &k);

}
