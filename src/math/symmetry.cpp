#include "symmetry.h"

#include <algorithm>
#include <cmath>
#include <iterator>
#include <set>
#include <stdexcept>
#include <tuple>

#include "../utils/error.h"

namespace librpa_int
{

bool same_fractional_kpoint(const Vector3_Order<double>& lhs,
                            const Vector3_Order<double>& rhs,
                            const double tol)
{
    return nearly_integer_vector(lhs - rhs, tol);
}

bool try_fold_fractional_kpoint_to_target(const Vector3_Order<double>& kpoint,
                                          const Vector3_Order<double>& target_kpoint,
                                          const double tol,
                                          FoldedKPoint& folded)
{
    const auto fold = kpoint - target_kpoint;
    if (!nearly_integer_vector(fold, tol))
    {
        return false;
    }
    folded.kpoint = target_kpoint;
    folded.fold_G = round_to_integer_vector(fold);
    return true;
}

static bool find_atom_to_target_mapping(const std::vector<Vector3_Order<double>>& atom_positions,
                                        const SpaceGroupSymOps& operations,
                                        const int atom,
                                        const int target,
                                        const double tol,
                                        AtomSymMapping& mapping)
{
    const auto atom_position = atom_positions[atom];
    const auto target_position = atom_positions[target];

    for (std::size_t isym = 0; isym < operations.size(); ++isym)
    {
        const auto& operation = operations[isym];
        const Vector3_Order<double> transformed =
            apply_space_group_symmetry_operation(operation, atom_position);
        const Vector3_Order<double> return_lattice = transformed - target_position;
        if (!nearly_integer_vector(return_lattice, tol))
        {
            continue;
        }
        mapping.inequivalent_atom = target;
        mapping.isym = static_cast<int>(isym);
        mapping.return_lattice = round_to_integer_vector(return_lattice);
        return true;
    }
    return false;
}

static bool find_kpoint_sym_mapping(const Vector3_Order<double>& representative_kpoint,
                                    const Vector3_Order<double>& member_kpoint,
                                    const SpaceGroupSymOps& operations,
                                    const double tol,
                                    KPointSymMapping& mapping)
{
    for (std::size_t isym = 0; isym < operations.size(); ++isym)
    {
        const auto rotated =
            apply_space_group_rotation_to_kpoint(operations[isym], representative_kpoint);
        FoldedKPoint folded;
        if (!try_fold_fractional_kpoint_to_target(rotated, member_kpoint, tol, folded))
        {
            continue;
        }
        mapping.isym = static_cast<int>(isym);
        mapping.fold_G = folded.fold_G;
        return true;
    }
    return false;
}

static int choose_kpoint_star_representative(
    const KPointStar& star,
    const std::size_t fallback_full_k_index,
    const std::vector<Vector3_Order<double>>& preferred_representative_kpoints,
    const double tol)
{
    for (const auto& preferred_kpoint : preferred_representative_kpoints)
    {
        for (std::size_t member_index = 0; member_index < star.members.size(); ++member_index)
        {
            FoldedKPoint folded;
            if (try_fold_fractional_kpoint_to_target(
                    preferred_kpoint, star.members[member_index].kpoint, tol, folded))
            {
                return static_cast<int>(member_index);
            }
        }
    }

    const auto representative_member = std::find_if(
        star.members.begin(),
        star.members.end(),
        [fallback_full_k_index](const KPointStarMember& member) {
            return member.full_k_index == static_cast<int>(fallback_full_k_index);
        });
    if (representative_member == star.members.end())
    {
        throw std::runtime_error("k-star does not contain its representative k-point");
    }
    return static_cast<int>(std::distance(star.members.begin(), representative_member));
}

bool preserves_lattice_metric(const Matrix3& rotation,
                              const Matrix3& lattice_vectors,
                              const double tol)
{
    // Test the metric in an orthonormal frame. An absolute comparison of
    // R G R^T and G has units of length squared: merely changing the lattice
    // unit or rounding a large cell can make a valid fractional operation
    // fail, causing callers to reinterpret it as a Cartesian rotation.
    const Matrix3 cartesian = row_fractional_rotation_to_cartesian(rotation, lattice_vectors);
    return is_same_matrix(cartesian.Transpose() * cartesian, Matrix3::IDENTITY, tol);
}

static Vector3_Order<int> rotate_rspace_vector(
    const Vector3_Order<int>& R,
    const SpaceGroupAtomMapping<int>& op_info,
    const SpaceGroupSymOp& op,
    const int atom_from_i,
    const int atom_from_j,
    const double tol)
{
    const auto to_double = [](const Vector3_Order<int>& vec) {
        return Vector3_Order<double>(static_cast<double>(vec.x),
                                     static_cast<double>(vec.y),
                                     static_cast<double>(vec.z));
    };
    const Vector3_Order<double> R_double{static_cast<double>(R.x),
                                         static_cast<double>(R.y),
                                         static_cast<double>(R.z)};
    const Vector3_Order<double> rotated_cell = multiply_row_vector(R_double, op.rotation);
    const Vector3_Order<double> common_shift = to_double(op_info.return_lattice[atom_from_i]);
    const Vector3_Order<double> mapped_j_cell =
        to_double(op_info.return_lattice[atom_from_j]) + rotated_cell - common_shift;
    if (!nearly_integer_vector(mapped_j_cell, tol))
    {
        throw std::runtime_error("Real-space symmetry generated a non-integer lattice vector");
    }
    return round_to_integer_vector(mapped_j_cell);
}

using RSpaceKey = std::tuple<int, int, Vector3_Order<int>>;

static bool rspace_representative_less(const RSpaceKey& lhs,
                                       const RSpaceKey& rhs)
{
    if (std::get<0>(lhs) != std::get<0>(rhs))
    {
        return std::get<0>(lhs) < std::get<0>(rhs);
    }
    if (std::get<1>(lhs) != std::get<1>(rhs))
    {
        return std::get<1>(lhs) < std::get<1>(rhs);
    }
    const auto lhs_norm = std::get<2>(lhs).norm_l1();
    const auto rhs_norm = std::get<2>(rhs).norm_l1();
    if (lhs_norm != rhs_norm)
    {
        return lhs_norm < rhs_norm;
    }
    return std::get<2>(lhs) < std::get<2>(rhs);
}

SpaceGroupSymOp compose_space_group_symmetry_operations(
    const SpaceGroupSymOp& lhs,
    const SpaceGroupSymOp& rhs)
{
    if (lhs.use_row_convention != rhs.use_row_convention)
    {
        throw std::invalid_argument("cannot compose symmetry operations with mixed conventions");
    }

    SpaceGroupSymOp composed;
    composed.use_row_convention = lhs.use_row_convention;
    if (composed.use_row_convention)
    {
        composed.rotation = multiply_space_group_rotation_matrices(lhs.rotation, rhs.rotation);
        composed.translation = multiply_row_vector(lhs.translation, rhs.rotation) + rhs.translation;
    }
    else
    {
        composed.rotation = multiply_space_group_rotation_matrices(rhs.rotation, lhs.rotation);
        composed.translation =
            Vector3_Order<double>(rhs.rotation * lhs.translation) + rhs.translation;
    }
    return composed;
}

Vector3_Order<double> apply_space_group_symmetry_operation(
    const SpaceGroupSymOp& operation,
    const Vector3_Order<double>& coord)
{
    if (operation.use_row_convention)
    {
        return multiply_row_vector(coord, operation.rotation) + operation.translation;
    }
    return Vector3_Order<double>(operation.rotation * coord) + operation.translation;
}

Vector3_Order<double> apply_space_group_rotation_to_kpoint(
    const SpaceGroupSymOp& operation,
    const Vector3_Order<double>& kpoint)
{
    const Matrix3 reciprocal_rotation = operation.rotation.Inverse().Transpose();
    if (operation.use_row_convention)
    {
        return multiply_row_vector(kpoint, reciprocal_rotation);
    }
    return Vector3_Order<double>(reciprocal_rotation * kpoint);
}

FoldedKPoint fold_fractional_kpoint_to_targets(
    const Vector3_Order<double>& kpoint,
    const std::vector<Vector3_Order<double>>& target_kpoints,
    const double tol)
{
    FoldedKPoint matched;
    bool found = false;
    for (std::size_t ik = 0; ik < target_kpoints.size(); ++ik)
    {
        FoldedKPoint candidate;
        if (!try_fold_fractional_kpoint_to_target(kpoint, target_kpoints[ik], tol, candidate))
        {
            continue;
        }
        if (found)
        {
            throw std::runtime_error("fractional k-point folds to more than one target k-point");
        }
        candidate.target_k_index = static_cast<int>(ik);
        matched = candidate;
        found = true;
    }
    if (!found)
    {
        throw std::runtime_error("fractional k-point does not fold to any target k-point");
    }
    return matched;
}

static std::vector<AtomSymMapping> build_fractional_atom_to_inequivalent_symmetry_mapping(
    const std::vector<Vector3_Order<double>>& atom_positions,
    const SpaceGroupSymOps& operations,
    const double tol)
{
    std::vector<AtomSymMapping> mappings(atom_positions.size());
    std::vector<int> inequivalent_atoms;

    for (std::size_t atom_index = 0; atom_index < atom_positions.size(); ++atom_index)
    {
        const int atom = static_cast<int>(atom_index);
        auto& mapping = mappings[atom_index];
        for (const int inequivalent_atom : inequivalent_atoms)
        {
            if (find_atom_to_target_mapping(atom_positions,
                                            operations,
                                            atom,
                                            inequivalent_atom,
                                            tol,
                                            mapping))
            {
                break;
            }
        }
        if (mapping.inequivalent_atom >= 0)
        {
            continue;
        }

        mapping.inequivalent_atom = atom;
        find_atom_to_target_mapping(atom_positions, operations, atom, atom, tol, mapping);
        inequivalent_atoms.push_back(atom);
    }

    return mappings;
}

const SpaceGroupSymOp SpaceGroupSymOp::IDENTITY{
    {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0}, {0.0, 0.0, 0.0}, true};

const SpaceGroupSymOp SpaceGroupSymOp::INVERSE{
    {-1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, -1.0}, {0.0, 0.0, 0.0}, true};

const SpaceGroupSymOp SpaceGroupSymOp::C41_Z{
    {0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0}, {0.0, 0.0, 0.0}, true};

int find_identity_symmetry_operation(
    const SpaceGroupSymOps &operations)
{
    for (std::size_t isym = 0; isym != operations.size(); ++isym)
    {
        const auto &op = operations.at(isym);
        if (op.is_identity())
        {
            return static_cast<int>(isym);
        }
    }
    throw LIBRPA_RUNTIME_ERROR("Symmetry operations do not contain the identity operation");
}

std::vector<AtomSymMapping> build_atom_to_inequivalent_symmetry_mapping(
    const std::vector<Vector3_Order<double>>& atom_positions_frac,
    const Matrix3& lattice_vectors,
    const SpaceGroupSymOps& fractional_operations,
    const double tol)
{
    if (std::fabs(lattice_vectors.Det()) < 1e-14)
    {
        throw std::invalid_argument("lattice vectors must be non-singular");
    }

    return build_fractional_atom_to_inequivalent_symmetry_mapping(
        atom_positions_frac, fractional_operations, tol);
}

std::vector<int> collect_inequivalent_atoms(
    const std::vector<AtomSymMapping>& mappings)
{
    std::vector<int> atoms;
    for (const auto& mapping : mappings)
    {
        if (mapping.inequivalent_atom < 0
            || std::find(atoms.begin(), atoms.end(), mapping.inequivalent_atom) != atoms.end())
        {
            continue;
        }
        atoms.push_back(mapping.inequivalent_atom);
    }
    return atoms;
}

SpaceGroupRSpaceSector build_space_group_rspace_irreducible_sector(
    const SpaceGroupSymOps& fractional_operations,
    const std::map<int, Vector3_Order<double>>& coord_frac,
    const std::map<int, int>& atom_to_type,
    const std::vector<Vector3_Order<int>>& Rlist,
    const Matrix3* lattice_vectors,
    const double atom_map_tol,
    const double coord_tol)
{
    if (fractional_operations.empty())
    {
        throw std::runtime_error("Real-space symmetry operations are unavailable");
    }
    if (atom_to_type.empty())
    {
        throw std::runtime_error("Atom-to-type mapping is unavailable for real-space symmetry");
    }

    std::vector<SpaceGroupAtomMapping<int>> op_infos;
    op_infos.reserve(fractional_operations.size());
    for (const auto& op : fractional_operations)
    {
        op_infos.push_back(get_space_group_atom_mapping(op, coord_frac, atom_to_type, atom_map_tol));
    }

    std::vector<bool> use_operation(fractional_operations.size(), true);
    if (lattice_vectors != nullptr)
    {
        for (std::size_t isym = 0; isym < fractional_operations.size(); ++isym)
        {
            use_operation[isym] =
                preserves_lattice_metric(fractional_operations[isym].rotation, *lattice_vectors);
        }
    }

    const std::set<Vector3_Order<int>> Rset(Rlist.begin(), Rlist.end());
    std::set<RSpaceKey> uncovered;
    for (int atom_i = 0; atom_i < static_cast<int>(atom_to_type.size()); ++atom_i)
    {
        for (int atom_j = 0; atom_j < static_cast<int>(atom_to_type.size()); ++atom_j)
        {
            for (const auto& R : Rlist)
            {
                uncovered.insert({atom_i, atom_j, R});
            }
        }
    }

    SpaceGroupRSpaceSector irreducible_sector;
    while (!uncovered.empty())
    {
        const auto seed = *uncovered.begin();
        const int seed_i = std::get<0>(seed);
        const int seed_j = std::get<1>(seed);
        const Vector3_Order<int> seed_R = std::get<2>(seed);

        std::set<RSpaceKey> orbit{seed};
        for (std::size_t isym = 0; isym < fractional_operations.size(); ++isym)
        {
            if (!use_operation[isym])
            {
                continue;
            }
            const auto& op_info = op_infos[isym];
            const auto full_I = op_info.atom_map[seed_i];
            const auto full_J = op_info.atom_map[seed_j];
            const auto full_R = rotate_rspace_vector(
                seed_R, op_info, fractional_operations[isym], seed_i, seed_j, coord_tol);
            if (Rset.count(full_R) != 0)
            {
                orbit.insert({full_I, full_J, full_R});
            }
        }

        const auto representative =
            *std::min_element(orbit.begin(), orbit.end(), rspace_representative_less);
        irreducible_sector[{std::get<0>(representative), std::get<1>(representative)}].insert(
            std::get<2>(representative));

        for (const auto& member : orbit)
        {
            uncovered.erase(member);
        }
    }

    return irreducible_sector;
}

std::vector<KPointStar> build_kpoint_stars(
    const std::vector<Vector3_Order<double>>& full_kpoints_frac,
    const SpaceGroupSymOps& fractional_operations,
    const double tol)
{
    return build_kpoint_stars(
        full_kpoints_frac, fractional_operations, std::vector<Vector3_Order<double>>{}, tol);
}

std::vector<KPointStar> build_kpoint_stars(
    const std::vector<Vector3_Order<double>>& full_kpoints_frac,
    const SpaceGroupSymOps& fractional_operations,
    const std::vector<Vector3_Order<double>>& preferred_representative_kpoints,
    const double tol)
{
    if (fractional_operations.empty())
    {
        throw std::invalid_argument("at least one symmetry operation is required to build k-stars");
    }

    const auto& kpoints = full_kpoints_frac;

    std::vector<KPointStar> stars;
    std::vector<bool> used(kpoints.size(), false);
    for (std::size_t ik_rep = 0; ik_rep < kpoints.size(); ++ik_rep)
    {
        if (used[ik_rep])
        {
            continue;
        }

        KPointStar star;
        const auto seed_kpoint = kpoints[ik_rep];

        // NOTE: O(n_k^2 * n_sym); add a grid hash if this becomes hot.
        for (std::size_t ik = 0; ik < kpoints.size(); ++ik)
        {
            if (used[ik])
            {
                continue;
            }
            for (std::size_t isym = 0; isym < fractional_operations.size(); ++isym)
            {
                const auto rotated =
                    apply_space_group_rotation_to_kpoint(fractional_operations[isym], seed_kpoint);
                FoldedKPoint folded;
                if (!try_fold_fractional_kpoint_to_target(rotated, kpoints[ik], tol, folded))
                {
                    continue;
                }

                KPointStarMember member;
                member.full_k_index = static_cast<int>(ik);
                member.kpoint = kpoints[ik];
                star.members.push_back(member);
                break;
            }
        }

        if (star.members.empty())
        {
            throw std::runtime_error("failed to build a k-star member for the representative k-point");
        }
        star.representative_k_index =
            choose_kpoint_star_representative(star, ik_rep, preferred_representative_kpoints, tol);
        const auto chosen_representative_kpoint =
            star.members[star.representative_k_index].kpoint;
        star.sym_mappings.reserve(star.members.size());
        for (const auto& member : star.members)
        {
            KPointSymMapping sym_mapping;
            if (!find_kpoint_sym_mapping(
                    chosen_representative_kpoint,
                    member.kpoint,
                    fractional_operations,
                    tol,
                    sym_mapping))
            {
                throw std::runtime_error("failed to build symmetry mapping for k-star member");
            }
            star.sym_mappings.push_back(sym_mapping);
        }
        for (const auto& member : star.members)
        {
            used[static_cast<std::size_t>(member.full_k_index)] = true;
        }
        stars.push_back(std::move(star));
    }

    return stars;
}

Matrix3 axial_rotation_of(const Matrix3 &cartesian_rotation)
{
    return cartesian_rotation * cartesian_rotation.Det();
}

std::array<std::complex<double>, 4> so3_to_su2(const Matrix3 &proper_rotation)
{
    using cplx = std::complex<double>;
    const Matrix3 &R = proper_rotation;
    const double trace = R.e11 + R.e22 + R.e33;

    // Antisymmetric part: v = 2 sin(theta) * n. With Matrix3 acting on column
    // vectors ((R u)_i = R_ij u_j), the axis is
    // n ~ (R_zy - R_yz, R_xz - R_zx, R_yx - R_xy) — e.g. the +90 deg rotation
    // about +z, R = [[0,-1,0],[1,0,0],[0,0,1]], gives n = +z and the physical
    // spinor U = diag(exp(-i pi/4), exp(i pi/4)), and the covering map
    // satisfies su2_to_so3(so3_to_su2(R)) == R. The planning report writes the
    // same tuple with swapped index order, (R_yz - R_zy, ...); that reading
    // corresponds to the row-vector action x' = x R and yields the adjoint
    // spinor U^dagger = U(-n, theta) for the column-vector matrix fed here.
    const double vx = R.e32 - R.e23;
    const double vy = R.e13 - R.e31;
    const double vz = R.e21 - R.e12;

    // cos(theta/2) = sqrt((Tr R + 1) / 4), non-negative for theta in [0, pi].
    // Near theta = pi this is ill-conditioned (d cos/d theta = 0): rounding
    // noise of the input shows up as sqrt-level noise here, so the pi branch
    // below must not trust it.
    const double half_cos = 0.5 * std::sqrt(std::max(trace + 1.0, 0.0));

    if (half_cos > 1e-7)
    {
        // Generic branch. The axis direction comes from the antisymmetric
        // part, n = v / |v|, and sin(theta/2) = sqrt(1 - cos(theta/2)^2) is
        // well-conditioned everywhere outside the pi neighborhood (unlike
        // |v| = 2 sin theta near theta = 0, where sin(theta/2) n = v/2 stays
        // smooth and the direction noise is suppressed by the tiny
        // sin(theta/2) factor). No 0/0 occurs at theta = 0.
        const double vnorm = std::sqrt(vx * vx + vy * vy + vz * vz);
        const double half_sin = std::sqrt(std::max(1.0 - half_cos * half_cos, 0.0));
        double sx = 0.0, sy = 0.0, sz = 0.0;
        if (vnorm > 1e-14)
        {
            const double scale = half_sin / vnorm;
            sx = vx * scale; // sin(theta/2) * n_x
            sy = vy * scale;
            sz = vz * scale;
        }
        // U = cos(theta/2) I - i sin(theta/2) (n . sigma)
        return {cplx(half_cos, -sz), cplx(-sy, -sx),
                cplx(sy, -sx), cplx(half_cos, sz)};
    }

    // theta = pi branch: R = 2 n n^T - I, so n_i^2 = (R_ii + 1) / 2.
    // Anchor the sign on the largest diagonal element and recover the other
    // components from R_ij + R_ji = 4 n_i n_j (ABACUS commit 1ed54ee23 fix).
    double nx = std::sqrt(std::max(0.5 * (R.e11 + 1.0), 0.0));
    double ny = std::sqrt(std::max(0.5 * (R.e22 + 1.0), 0.0));
    double nz = std::sqrt(std::max(0.5 * (R.e33 + 1.0), 0.0));
    if (nx >= ny && nx >= nz)
    {
        ny = (R.e12 + R.e21) / (4.0 * nx);
        nz = (R.e13 + R.e31) / (4.0 * nx);
    }
    else if (ny >= nz)
    {
        nx = (R.e12 + R.e21) / (4.0 * ny);
        nz = (R.e23 + R.e32) / (4.0 * ny);
    }
    else
    {
        nx = (R.e13 + R.e31) / (4.0 * nz);
        ny = (R.e23 + R.e32) / (4.0 * nz);
    }
    // U = -i (n . sigma); the real part cos(theta/2) is set to zero because
    // its trace-based estimate is ill-conditioned in this branch (see above).
    return {cplx(0.0, -nz), cplx(-ny, -nx),
            cplx(ny, -nx), cplx(0.0, nz)};
}

namespace
{

using su2_matrix_t = std::array<std::complex<double>, 4>;

su2_matrix_t su2_multiply(const su2_matrix_t &a, const su2_matrix_t &b)
{
    return {a[0] * b[0] + a[1] * b[2], a[0] * b[1] + a[1] * b[3],
            a[2] * b[0] + a[3] * b[2], a[2] * b[1] + a[3] * b[3]};
}

su2_matrix_t su2_dagger(const su2_matrix_t &a)
{
    return {std::conj(a[0]), std::conj(a[2]), std::conj(a[1]), std::conj(a[3])};
}

std::complex<double> su2_trace(const su2_matrix_t &a)
{
    return a[0] + a[3];
}

} // namespace

Matrix3 su2_to_so3(const std::array<std::complex<double>, 4> &U)
{
    const std::complex<double> imag_i(0.0, 1.0);
    const su2_matrix_t sigma[3] = {
        {0.0, 1.0, 1.0, 0.0},
        {0.0, -imag_i, imag_i, 0.0},
        {1.0, 0.0, 0.0, -1.0}};
    const su2_matrix_t U_dagger = su2_dagger(U);
    std::array<double, 9> w{};
    for (int i = 0; i != 3; ++i)
    {
        for (int j = 0; j != 3; ++j)
        {
            w[3 * i + j] = 0.5 * std::real(
                su2_trace(su2_multiply(su2_multiply(sigma[i], U),
                                       su2_multiply(sigma[j], U_dagger))));
        }
    }
    return Matrix3(w);
}

} // namespace librpa_int
