/*!
 * @file symmetry_types.h
 * @brief Shared symmetry data structures.
 */
#pragma once

#include <array>
#include <complex>
#include <cstddef>
#include <map>
#include <set>
#include <vector>

#include "atom.h"
#include "../math/complexmatrix.h"
#include "../math/symmetry.h"
#include "../math/vector3_order.h"

namespace librpa_int
{

using symmetry_R_t = std::array<int, 3>;
using symmetry_irreducible_sector_t = std::map<atpair_t, std::set<symmetry_R_t>>;

using SymmetryOperation = SpaceGroupSymOp;

/*!
 * @brief Physical origin of the SU(2) spin action attached to a spatial operation.
 */
enum class SymmetrySpinActionSource
{
    //! Ordinary space group / grey group: U_s is the identity.
    Identity,
    //! U_s is explicitly provided by the upstream spin-space-group input.
    ExplicitSpinSpace,
    //! U_s = U[det(Q) Q] reconstructed from the spatial rotation (SOC magnetic group).
    DerivedFromSpatialSOC
};

/*!
 * @brief Storage layout of the spin degree of freedom.
 *
 * Declared in Phase 1 as metadata only; not wired to any production path yet.
 */
enum class SpinStorageMode
{
    Scalar,
    CollinearChannels,
    Spinor2
};

/*!
 * @brief Physical symmetry mode of the calculation.
 *
 * Declared in Phase 1 as metadata only; not wired to any production path yet.
 */
enum class SymmetryPhysicsMode
{
    SpatialOnly,
    SpinSpaceGroup,
    MagneticSpaceGroupSOC
};

/*!
 * @brief Classification of the pure-spin stabilizer of a spin texture.
 *
 * Declared in Phase 1 as metadata only; not wired to any production path yet.
 */
enum class PureSpinStabilizerKind
{
    None,
    SU2Invariant,
    U1Axis,
    CoplanarZ2,
    ExplicitFinite
};

/*!
 * @brief Pure-spin stabilizer metadata: rotations that act only on spin.
 *
 * Declared in Phase 1 as metadata only; generators stay empty until pure-spin
 * operations are accepted as input.
 */
struct PureSpinStabilizer
{
    PureSpinStabilizerKind kind = PureSpinStabilizerKind::None;
    Vector3_Order<double> axis_or_normal{0.0, 0.0, 1.0};
    std::vector<std::array<std::complex<double>, 4>> generators;
};

/*!
 * @brief One full spin-space operation (g, U_s, eta).
 *
 * `spatial_id` indexes the spatial operation pool, `spin_u` is the SU(2) spin
 * rotation stored as a row-major 2x2 matrix, and `antiunitary` flags the
 * antiunitary (time-reversal-like) part. For ordinary space-group input the
 * spin action is the identity and the source is Identity.
 */
struct SymmetrySpinOperation
{
    std::size_t spatial_id = 0;
    std::array<std::complex<double>, 4> spin_u{1.0, 0.0, 0.0, 1.0};
    bool antiunitary = false;
    SymmetrySpinActionSource spin_source = SymmetrySpinActionSource::Identity;
    //! Preserve an upstream SOC matrix instead of deriving it from the lattice.
    bool spin_u_provided = false;
};

/*!
 * @brief One deduplicated geometric action on k/q points.
 *
 * For k/q-space the geometry key is (spatial_id, antiunitary); the canonical
 * operation is the first spin operation with that key and
 * `equivalent_operation_ids` collects all spin operations sharing the geometry
 * (pure-spin duplicates included, so they never widen a star).
 */
struct SymmetryGeometricAction
{
    std::size_t canonical_operation_id = 0;
    std::vector<std::size_t> equivalent_operation_ids;
};

/*!
 * @brief Cached metadata of one spatial operation in the operation pool.
 *
 * `reciprocal_rotation` is the dual rotation applied to fractional k/q points,
 * identical to the one used by apply_space_group_rotation_to_kpoint.
 * `atom_map` and `return_lattice` are lazy cache slots filled on demand.
 */
struct SymmetrySpatialOperationRecord
{
    SpaceGroupSymOp spatial;
    Matrix3 reciprocal_rotation;
    std::vector<atom_t> atom_map;
    std::vector<Vector3_Order<int>> return_lattice;
};

/*!
 * @brief Atom-resolved k-space symmetry information exported by a symmetry convention.
 */
struct SymmetryKAtomRotation
{
    int atom_from = -1;
    int atom_to = -1;
    int atom_type = -1;
    int lmax = -1;
    std::map<int, ComplexMatrix> bloch_rsh_rotations;
};

/*!
 * @brief One member of an irreducible k-star.
 */
struct SymmetryKStarMember
{
    int spatial_isym = -1;
    bool time_reversal = false;
    Vector3_Order<double> k_bz{0.0, 0.0, 0.0};
    //! Index into SymmetryContext::kspace_actions; derived from
    //! (spatial_isym, time_reversal) and kept consistent with them.
    std::size_t action_id = 0;
    std::vector<SymmetryKAtomRotation> atom_rotations;
};

/*!
 * @brief One irreducible k-star exported by a symmetry convention.
 */
struct SymmetryKStar
{
    int star_index = -1;
    Vector3_Order<double> k_ibz{0.0, 0.0, 0.0};
    std::vector<SymmetryKStarMember> members;
};

/*!
 * @brief Explicit mapping between one loaded LibRPA IBZ q-index and one symmetry k-star.
 *
 * The mapping stores the full-BZ q keys that LibRPA should use for every star member.
 * Symmetry coordinates are treated as the source of truth; when LibRPA already
 * has an equivalent internal q key, that exact key is reused to keep later lookups
 * aligned with existing storage.
 */
struct SymmetryKStarGridMappingEntry
{
    int iq_ibz = -1;
    int star_list_index = -1;
    std::vector<Vector3_Order<double>> member_q_bz_keys;
};

/*!
 * @brief One full-BZ k-point member expanded from a symmetry IBZ k-star.
 */
struct SymmetryFullKpointMemberEntry
{
    int ik_ibz = -1;
    int star_list_index = -1;
    int member_index = -1;
    int spatial_isym = -1;
    bool time_reversal = false;
    Vector3_Order<double> k_bz{0.0, 0.0, 0.0};
    //! Index into SymmetryContext::kspace_actions, copied from the star member.
    std::size_t action_id = 0;
};

/*!
 * @brief One full real-space member generated from an irreducible {atom pair, R}.
 *
 * `isym` indexes `SymmetryContext::rspace_operations` (the spatial part) and is
 * kept for all pre-magnetic consumers. `operation_id` indexes
 * `SymmetryContext::spin_operations` and records which full (g, U_s, eta)
 * operation generated this member; it equals `kSymmetryRSpaceOperationIdNone`
 * when the star was built without spin-operation metadata (legacy contexts).
 */
struct SymmetryRSpaceRestoreMember
{
    static constexpr std::size_t kOperationIdNone = ~static_cast<std::size_t>(0);
    int isym = -1;
    atpair_t full_atom_pair;
    Vector3_Order<int> full_R{0, 0, 0};
    std::size_t operation_id = kOperationIdNone;
};

using symmetry_rspace_sector_stars_t =
    std::map<atpair_t, std::map<Vector3_Order<int>, std::vector<SymmetryRSpaceRestoreMember>>>;
using symmetry_atom_block_matrix_map_t = std::map<atom_t, std::map<atom_t, ComplexMatrix>>;
//! One spin-channel real-space tensor map keyed {I, {J, R}} with dense AO blocks.
using symmetry_rspace_block_map_t =
    std::map<int, std::map<std::pair<int, std::array<int, 3>>, ComplexMatrix>>;
using symmetry_kstar_member_kfrac_targets_t =
    std::vector<std::vector<Vector3_Order<double>>>;
using symmetry_kstar_representative_indices_t = std::vector<int>;

} // namespace librpa_int
