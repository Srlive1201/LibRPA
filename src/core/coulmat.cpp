#include <algorithm>
#include <cmath>
#include <complex>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include "../utils/constants.h"
#include "../io/global_io.h"
#include "../math/utils_matrix_mpi.h"
#include "symmetry_context.h"
#include "coulmat.h"

namespace librpa_int {

static bool are_equivalent_symmetry_qpoints(const Vector3_Order<double>& lhs,
                                            const Vector3_Order<double>& rhs,
                                            const double tol = 1e-5)
{
    const auto same_component = [tol](const double lhs_component, const double rhs_component) {
        return std::abs((lhs_component - rhs_component) - std::round(lhs_component - rhs_component))
               < tol;
    };
    return same_component(lhs.x, rhs.x) && same_component(lhs.y, rhs.y)
           && same_component(lhs.z, rhs.z);
}

template <typename QMap>
static typename QMap::const_iterator find_matching_symmetry_qpoint(
    const QMap& q_map,
    const Vector3_Order<double>& q_target)
{
    const auto exact_iter = q_map.find(q_target);
    if (exact_iter != q_map.end())
    {
        return exact_iter;
    }

    return std::find_if(q_map.begin(), q_map.end(), [&q_target](const auto& entry) {
        return are_equivalent_symmetry_qpoints(entry.first, q_target);
    });
}

static std::map<atom_t, size_t> build_atom_nabf_map(const AtomicBasis& basis_abf)
{
    std::map<atom_t, size_t> atom_nabf;
    for (atom_t atom = 0; atom != static_cast<atom_t>(basis_abf.n_atoms); ++atom)
    {
        atom_nabf[atom] = basis_abf[static_cast<int>(atom)];
    }
    return atom_nabf;
}

static bool has_complete_symmetry_abf_ibz_coverage(
    const atpair_k_cplx_mat_t& blocks_by_q_ibz,
    const std::map<atom_t, size_t>& atom_nabf,
    const PeriodicBoundaryData& pbc)
{
    for (std::size_t atom_i = 0; atom_i < atom_nabf.size(); ++atom_i)
    {
        for (std::size_t atom_j = atom_i; atom_j < atom_nabf.size(); ++atom_j)
        {
            const auto upper_iter = blocks_by_q_ibz.find(static_cast<atom_t>(atom_i));
            const bool has_upper =
                upper_iter != blocks_by_q_ibz.end()
                && upper_iter->second.count(static_cast<atom_t>(atom_j)) != 0;

            const auto lower_iter = blocks_by_q_ibz.find(static_cast<atom_t>(atom_j));
            const bool has_lower =
                lower_iter != blocks_by_q_ibz.end()
                && lower_iter->second.count(static_cast<atom_t>(atom_i)) != 0;

            if (!has_upper && !has_lower)
            {
                return false;
            }

            const auto& q_blocks =
                has_upper ? upper_iter->second.at(static_cast<atom_t>(atom_j))
                          : lower_iter->second.at(static_cast<atom_t>(atom_i));
            for (const auto& q_ibz : pbc.klist)
            {
                if (find_matching_symmetry_qpoint(q_blocks, q_ibz) == q_blocks.end())
                {
                    return false;
                }
            }
        }
    }
    return true;
}

static librpa_int::symmetry_atom_block_matrix_map_t collect_symmetry_abf_ibz_blocks_for_q(
    const atpair_k_cplx_mat_t& blocks_by_q,
    const Vector3_Order<double>& q_ibz_internal)
{
    librpa_int::symmetry_atom_block_matrix_map_t blocks_ibz;
    for (const auto& atom_i_pair : blocks_by_q)
    {
        const auto atom_i = atom_i_pair.first;
        for (const auto& atom_j_pair : atom_i_pair.second)
        {
            const auto atom_j = atom_j_pair.first;
            const auto q_iter = find_matching_symmetry_qpoint(atom_j_pair.second, q_ibz_internal);
            if (q_iter != atom_j_pair.second.end())
            {
                blocks_ibz[atom_i][atom_j] = *q_iter->second;
            }
        }
    }
    return blocks_ibz;
}

static std::vector<int> build_symmetry_atom_offsets(const std::map<atom_t, size_t>& atom_nabf)
{
    std::vector<int> offsets(atom_nabf.size() + 1, 0);
    for (std::size_t atom = 0; atom < atom_nabf.size(); ++atom)
    {
        offsets[atom + 1] = offsets[atom]
                            + static_cast<int>(atom_nabf.at(static_cast<atom_t>(atom)));
    }
    return offsets;
}

static ComplexMatrix build_dense_symmetry_hermitian_matrix_from_local_blocks(
    const librpa_int::symmetry_atom_block_matrix_map_t& local_blocks,
    const std::map<atom_t, size_t>& atom_nabf)
{
    const auto offsets = build_symmetry_atom_offsets(atom_nabf);
    ComplexMatrix dense(offsets.back(), offsets.back());
    for (const auto& atom_i_pair : local_blocks)
    {
        const int row_offset = offsets.at(static_cast<std::size_t>(atom_i_pair.first));
        const int expected_nrows = static_cast<int>(atom_nabf.at(atom_i_pair.first));
        for (const auto& atom_j_pair : atom_i_pair.second)
        {
            const int col_offset = offsets.at(static_cast<std::size_t>(atom_j_pair.first));
            const int expected_ncols = static_cast<int>(atom_nabf.at(atom_j_pair.first));
            const auto& block = atom_j_pair.second;
            if (block.nr != expected_nrows || block.nc != expected_ncols)
            {
                std::ostringstream oss;
                oss << "Dense V(q) symmetry restore block dimension mismatch for atom pair ("
                    << atom_i_pair.first << "," << atom_j_pair.first << "): block=" << block.nr
                    << "x" << block.nc << ", expected=" << expected_nrows << "x"
                    << expected_ncols;
                throw std::runtime_error(oss.str());
            }
            for (int row = 0; row < block.nr; ++row)
            {
                for (int col = 0; col < block.nc; ++col)
                {
                    const auto value = block(row, col);
                    dense(row_offset + row, col_offset + col) = value;
                    if (atom_i_pair.first != atom_j_pair.first)
                    {
                        dense(col_offset + col, row_offset + row) = std::conj(value);
                    }
                }
            }
        }
    }
    return dense;
}

static librpa_int::symmetry_atom_block_matrix_map_t build_symmetry_blocks_from_dense_matrix(
    const ComplexMatrix& dense_matrix,
    const std::map<atom_t, size_t>& atom_nabf)
{
    const auto offsets = build_symmetry_atom_offsets(atom_nabf);
    librpa_int::symmetry_atom_block_matrix_map_t atom_blocks;
    for (std::size_t atom_i = 0; atom_i < atom_nabf.size(); ++atom_i)
    {
        const int row_offset = offsets.at(atom_i);
        const int nrows = static_cast<int>(atom_nabf.at(static_cast<atom_t>(atom_i)));
        for (std::size_t atom_j = atom_i; atom_j < atom_nabf.size(); ++atom_j)
        {
            const int col_offset = offsets.at(atom_j);
            const int ncols = static_cast<int>(atom_nabf.at(static_cast<atom_t>(atom_j)));
            ComplexMatrix block(nrows, ncols);
            for (int row = 0; row < nrows; ++row)
            {
                for (int col = 0; col < ncols; ++col)
                {
                    block(row, col) = dense_matrix(row_offset + row, col_offset + col);
                }
            }
            atom_blocks[static_cast<atom_t>(atom_i)][static_cast<atom_t>(atom_j)] =
                std::move(block);
        }
    }
    return atom_blocks;
}

static librpa_int::symmetry_atom_block_matrix_map_t gather_symmetry_ibz_blocks_collective(
    const MpiCommHandler& comm_h,
    const librpa_int::symmetry_atom_block_matrix_map_t& blocks_ibz_local,
    const std::map<atom_t, size_t>& atom_nabf)
{
    if (comm_h.nprocs <= 1)
    {
        return blocks_ibz_local;
    }

    const auto dense_local =
        build_dense_symmetry_hermitian_matrix_from_local_blocks(blocks_ibz_local, atom_nabf);
    ComplexMatrix dense_global(dense_local.nr, dense_local.nc);
    allreduce_ComplexMatrix(dense_local, dense_global, comm_h.comm);
    return build_symmetry_blocks_from_dense_matrix(dense_global, atom_nabf);
}

static librpa_int::symmetry_irreducible_sector_t filter_symmetry_irreducible_sector_by_rlist(
    const librpa_int::symmetry_irreducible_sector_t& irreducible_sector,
    const std::vector<Vector3_Order<int>>& Rlist)
{
    librpa_int::symmetry_irreducible_sector_t filtered_sector;
    const std::set<Vector3_Order<int>> requested_rset(Rlist.begin(), Rlist.end());
    for (const auto& pair_Rs : irreducible_sector)
    {
        for (const auto& R_array : pair_Rs.second)
        {
            const Vector3_Order<int> R{R_array[0], R_array[1], R_array[2]};
            if (requested_rset.count(R) == 0)
            {
                continue;
            }
            filtered_sector[pair_Rs.first].insert(R_array);
        }
    }
    return filtered_sector;
}

static std::set<std::pair<atom_t, atom_t>> build_symmetry_irreducible_target_atom_pairs(
    const librpa_int::symmetry_irreducible_sector_t& irreducible_sector)
{
    std::set<std::pair<atom_t, atom_t>> target_atom_pairs;
    for (const auto& pair_Rs : irreducible_sector)
    {
        if (!pair_Rs.second.empty())
        {
            target_atom_pairs.insert(pair_Rs.first);
        }
    }
    return target_atom_pairs;
}

static std::set<std::pair<atom_t, atom_t>> collect_symmetry_local_target_atom_pairs(
    const atpair_k_cplx_mat_t& blocks_by_q,
    const bool return_ordered_atom_pair)
{
    std::set<std::pair<atom_t, atom_t>> target_atom_pairs;
    for (const auto& atom_i_pair : blocks_by_q)
    {
        for (const auto& atom_j_pair : atom_i_pair.second)
        {
            if (atom_j_pair.second.empty())
            {
                continue;
            }
            target_atom_pairs.insert({atom_i_pair.first, atom_j_pair.first});
            if (return_ordered_atom_pair)
            {
                target_atom_pairs.insert({atom_j_pair.first, atom_i_pair.first});
            }
        }
    }
    return target_atom_pairs;
}

static std::complex<double> build_ft_vq_phase(const PeriodicBoundaryData& pbc,
                                              const Vector3_Order<double>& q_internal,
                                              const Vector3_Order<int>& R)
{
    const auto q_frac = pbc.latvec * q_internal;
    const double ang = -(q_frac * R) * TWO_PI;
    return std::complex<double>(std::cos(ang), std::sin(ang))
           / static_cast<double>(pbc.get_n_cells_bvk());
}

static atpair_R_mat_t accumulate_symmetry_abf_irreducible_sector_vr(
    const MpiCommHandler& comm_h,
    const SymmetryContext& ctx,
    const AtomicBasis& basis_abf,
    const atpair_k_cplx_mat_t& blocks_by_q_ibz,
    const PeriodicBoundaryData& pbc)
{
    atpair_R_mat_t blocks_by_R_real;
    atpair_R_cplx_mat_t blocks_by_R_complex;
    const auto atom_nabf = build_atom_nabf_map(basis_abf);
    const auto abf_layouts = basis_abf.build_species_basis_layouts(ctx.atom_to_type);
    if (!symmetry_species_layouts_match_atom_counts(abf_layouts, ctx.atom_to_type, atom_nabf))
    {
        throw std::runtime_error("Auxiliary basis shell layout is inconsistent with atom_nabf");
    }
    const auto filtered_sector =
        filter_symmetry_irreducible_sector_by_rlist(ctx.irreducible_sector, pbc.Rlist);
    if (filtered_sector.empty())
    {
        return blocks_by_R_real;
    }

    const auto target_atom_pairs = build_symmetry_irreducible_target_atom_pairs(filtered_sector);
    for (const auto& pair_Rs : filtered_sector)
    {
        const auto atom_i = pair_Rs.first.first;
        const auto atom_j = pair_Rs.first.second;
        const int n_i = static_cast<int>(atom_nabf.at(atom_i));
        const int n_j = static_cast<int>(atom_nabf.at(atom_j));
        for (const auto& R_array : pair_Rs.second)
        {
            const Vector3_Order<int> R{R_array[0], R_array[1], R_array[2]};
            blocks_by_R_complex[atom_i][atom_j][R] =
                std::make_shared<ComplexMatrix>(n_i, n_j);
        }
    }

    for (const auto& star_mapping : ctx.kstar_grid_mapping)
    {
        const auto& star = ctx.kstars.at(static_cast<std::size_t>(star_mapping.star_list_index));

        const auto q_ibz_internal = pbc.klist.at(static_cast<std::size_t>(star_mapping.iq_ibz));
        auto blocks_ibz =
            collect_symmetry_abf_ibz_blocks_for_q(blocks_by_q_ibz, q_ibz_internal);
        blocks_ibz = gather_symmetry_ibz_blocks_collective(comm_h, blocks_ibz, atom_nabf);
        if (blocks_ibz.empty())
        {
            continue;
        }
        if (star.members.size() != star_mapping.member_q_bz_keys.size())
        {
            throw std::runtime_error(
                "Symmetry q-star mapping is inconsistent with the loaded full-q keys");
        }

        for (std::size_t imember = 0; imember < star.members.size(); ++imember)
        {
            const auto& member = star.members[imember];
            const auto q_bz_target_frac_vec =
                pbc.latvec * star_mapping.member_q_bz_keys[imember];
            const Vector3_Order<double> q_bz_target_frac{
                q_bz_target_frac_vec.x, q_bz_target_frac_vec.y, q_bz_target_frac_vec.z};
            librpa_int::symmetry_atom_block_matrix_map_t rotated_blocks;
            try
            {
                rotated_blocks = librpa_int::rotate_symmetry_kspace_operator_blocks(
                    ctx, abf_layouts, member, blocks_ibz, atom_nabf, star.k_ibz,
                    member.time_reversal, &target_atom_pairs, &q_bz_target_frac);
            }
            catch (const std::exception& ex)
            {
                std::ostringstream oss;
                oss << "Symmetry irreducible-sector FT failed for star=" << star.star_index
                    << ", member=" << imember << ", spatial_isym=" << member.spatial_isym
                    << ", time_reversal=" << (member.time_reversal ? "true" : "false") << ": "
                    << ex.what();
                throw std::runtime_error(oss.str());
            }

            const auto& q_internal = star_mapping.member_q_bz_keys[imember];
            for (const auto& atom_i_pair : rotated_blocks)
            {
                for (const auto& atom_j_pair : atom_i_pair.second)
                {
                    const auto sector_iter =
                        filtered_sector.find({atom_i_pair.first, atom_j_pair.first});
                    if (sector_iter == filtered_sector.end())
                    {
                        continue;
                    }

                    for (const auto& R_array : sector_iter->second)
                    {
                        const Vector3_Order<int> R{R_array[0], R_array[1], R_array[2]};
                        const auto phase = build_ft_vq_phase(pbc, q_internal, R);
                        *blocks_by_R_complex.at(atom_i_pair.first).at(atom_j_pair.first).at(R) +=
                            atom_j_pair.second * phase;
                    }
                }
            }
        }
    }

    for (const auto& atom_i_pair : blocks_by_R_complex)
    {
        for (const auto& atom_j_pair : atom_i_pair.second)
        {
            for (const auto& R_block : atom_j_pair.second)
            {
                blocks_by_R_real[atom_i_pair.first][atom_j_pair.first][R_block.first] =
                    std::make_shared<matrix>(R_block.second->real());
            }
        }
    }

    return blocks_by_R_real;
}

static atpair_R_mat_t accumulate_symmetry_abf_full_sector_vr(
    const MpiCommHandler& comm_h,
    const SymmetryContext& ctx,
    const AtomicBasis& basis_abf,
    const atpair_k_cplx_mat_t& blocks_by_q_ibz,
    const PeriodicBoundaryData& pbc,
    const bool return_ordered_atom_pair)
{
    atpair_R_mat_t blocks_by_R_real;
    atpair_R_cplx_mat_t blocks_by_R_complex;
    const auto atom_nabf = build_atom_nabf_map(basis_abf);
    const auto abf_layouts = basis_abf.build_species_basis_layouts(ctx.atom_to_type);
    if (!symmetry_species_layouts_match_atom_counts(abf_layouts, ctx.atom_to_type, atom_nabf))
    {
        throw std::runtime_error("Auxiliary basis shell layout is inconsistent with atom_nabf");
    }

    const auto target_atom_pairs =
        collect_symmetry_local_target_atom_pairs(blocks_by_q_ibz, return_ordered_atom_pair);
    for (const auto& atom_pair : target_atom_pairs)
    {
        const int n_i = static_cast<int>(atom_nabf.at(atom_pair.first));
        const int n_j = static_cast<int>(atom_nabf.at(atom_pair.second));
        for (const auto& R : pbc.Rlist)
        {
            blocks_by_R_complex[atom_pair.first][atom_pair.second][R] =
                std::make_shared<ComplexMatrix>(n_i, n_j);
        }
    }

    const auto full_grid_member_targets =
        build_symmetry_full_grid_kstar_member_kfrac_targets(ctx, pbc.kfrac_list);
    const bool use_full_grid_member_targets =
        full_grid_member_targets.size() == ctx.kstars.size();

    for (const auto& star_mapping : ctx.kstar_grid_mapping)
    {
        const auto& star = ctx.kstars.at(static_cast<std::size_t>(star_mapping.star_list_index));
        const auto q_ibz_internal = pbc.klist.at(static_cast<std::size_t>(star_mapping.iq_ibz));
        const auto q_ibz_frac = pbc.kfrac_list.at(static_cast<std::size_t>(star_mapping.iq_ibz));
        auto blocks_ibz =
            collect_symmetry_abf_ibz_blocks_for_q(blocks_by_q_ibz, q_ibz_internal);
        blocks_ibz = gather_symmetry_ibz_blocks_collective(comm_h, blocks_ibz, atom_nabf);
        if (target_atom_pairs.empty() || blocks_ibz.empty())
        {
            continue;
        }

        const auto rotation_atom_pairs =
            build_symmetry_upper_atom_pair_closure(star, target_atom_pairs);
        blocks_ibz = symmetrize_symmetry_ibz_kspace_operator_blocks(
            ctx, abf_layouts, q_ibz_frac, blocks_ibz, atom_nabf, &rotation_atom_pairs);
        if (star.members.size() != star_mapping.member_q_bz_keys.size())
        {
            throw std::runtime_error(
                "Symmetry q-star mapping is inconsistent with the loaded full-q keys");
        }

        for (std::size_t imember = 0; imember < star.members.size(); ++imember)
        {
            const auto& member = star.members[imember];
            const Vector3_Order<double> raw_q_bz_target_frac =
                use_full_grid_member_targets
                    ? full_grid_member_targets.at(
                          static_cast<std::size_t>(star_mapping.star_list_index)).at(imember)
                    : Vector3_Order<double>{pbc.latvec *
                                            star_mapping.member_q_bz_keys.at(imember)};
            const auto q_bz_target_frac =
                restrict_fractional_coordinate(raw_q_bz_target_frac);

            symmetry_atom_block_matrix_map_t rotated_blocks;
            try
            {
                rotated_blocks = rotate_symmetry_kspace_operator_blocks(
                    ctx, abf_layouts, member, blocks_ibz, atom_nabf, star.k_ibz,
                    member.time_reversal, &target_atom_pairs, &q_bz_target_frac);
            }
            catch (const std::exception& ex)
            {
                std::ostringstream oss;
                oss << "Symmetry full-sector V(q)->V(R) accumulation failed for star="
                    << star.star_index << ", member=" << imember
                    << ", spatial_isym=" << member.spatial_isym
                    << ", time_reversal=" << (member.time_reversal ? "true" : "false")
                    << ": " << ex.what();
                throw std::runtime_error(oss.str());
            }

            const auto& q_internal = star_mapping.member_q_bz_keys.at(imember);
            for (const auto& atom_i_pair : rotated_blocks)
            {
                for (const auto& atom_j_pair : atom_i_pair.second)
                {
                    for (const auto& R : pbc.Rlist)
                    {
                        const auto phase = build_ft_vq_phase(pbc, q_internal, R);
                        *blocks_by_R_complex.at(atom_i_pair.first)
                             .at(atom_j_pair.first)
                             .at(R) += atom_j_pair.second * phase;
                    }
                }
            }
        }
    }

    for (const auto& atom_i_pair : blocks_by_R_complex)
    {
        for (const auto& atom_j_pair : atom_i_pair.second)
        {
            for (const auto& R_block : atom_j_pair.second)
            {
                blocks_by_R_real[atom_i_pair.first][atom_j_pair.first][R_block.first] =
                    std::make_shared<matrix>(R_block.second->real());
            }
        }
    }
    return blocks_by_R_real;
}

static bool can_use_symmetry_qstar_ft_vq(const SymmetryContext& ctx,
                                         const AtomicBasis& basis_abf,
                                         const atpair_k_cplx_mat_t& coulmat_k,
                                         const PeriodicBoundaryData& pbc,
                                         const MpiCommHandler& comm_h)
{
    const auto atom_nabf = build_atom_nabf_map(basis_abf);
    const auto abf_layouts = basis_abf.build_species_basis_layouts(ctx.atom_to_type);
    if (!symmetry_species_layouts_match_atom_counts(abf_layouts, ctx.atom_to_type, atom_nabf))
    {
        return false;
    }
    const bool has_complete_or_distributed_coverage =
        comm_h.nprocs > 1
        || has_complete_symmetry_abf_ibz_coverage(coulmat_k, atom_nabf, pbc);
    return ctx.available
           && !ctx.kstars.empty()
           && ctx.kstars.size() == pbc.kfrac_list.size()
           && ctx.kstar_grid_mapping.size() == ctx.kstars.size()
           && ctx.count_kstar_members() == static_cast<std::size_t>(pbc.get_n_cells_bvk())
           && atom_nabf.size() == ctx.atom_to_type.size()
           && ctx.input_coord_frac.size() == atom_nabf.size()
           && pbc.klist.size() < static_cast<std::size_t>(pbc.get_n_cells_bvk())
           && has_complete_or_distributed_coverage
           && !ctx.rspace_operations.empty();
}

static bool can_use_symmetry_irreducible_sector_ft_vq(const SymmetryContext& ctx,
                                                      const AtomicBasis& basis_abf,
                                                      const atpair_k_cplx_mat_t& coulmat_k,
                                                      const PeriodicBoundaryData& pbc,
                                                      const MpiCommHandler& comm_h)
{
    return can_use_symmetry_qstar_ft_vq(ctx, basis_abf, coulmat_k, pbc, comm_h)
           && !ctx.irreducible_sector.empty()
           && !ctx.rspace_operations.empty();
}

atpair_R_mat_t FT_Vq(const MpiCommHandler &comm_h,
                     const AtomicBasis &basis_abf,
                     const SymmetryContext &symmetry_context,
                     const atpair_k_cplx_mat_t &coulmat_k,
                     const PeriodicBoundaryData &pbc,
                     bool return_ordered_atom_pair,
                     const bool use_symmetry_context,
                     bool* result_is_replicated)
{
    atpair_R_mat_t coulmat_R;
    if (result_is_replicated) *result_is_replicated = false;

    const auto &Rlist = pbc.Rlist;
    const auto &latvec = pbc.latvec;
    const auto &map_irk_ks = pbc.map_irk_ks;
    const auto n_k_points = pbc.get_n_cells_bvk();

    if (use_symmetry_context
        && can_use_symmetry_irreducible_sector_ft_vq(
            symmetry_context, basis_abf, coulmat_k, pbc, comm_h))
    {
        global::lib_printf_root(
            "EXX symmetry accumulates irreducible-sector `V(R)` directly from IBZ q-stars\n");
        if (result_is_replicated) *result_is_replicated = true;
        return accumulate_symmetry_abf_irreducible_sector_vr(
            comm_h, symmetry_context, basis_abf, coulmat_k, pbc);
    }

    if (use_symmetry_context
        && can_use_symmetry_qstar_ft_vq(
            symmetry_context, basis_abf, coulmat_k, pbc, comm_h))
    {
        global::lib_printf_root(
            "EXX symmetry accumulates full-sector `V(R)` directly from IBZ q-stars\n");
        return accumulate_symmetry_abf_full_sector_vr(
            comm_h, symmetry_context, basis_abf, coulmat_k, pbc,
            return_ordered_atom_pair);
    }

    for (auto R: Rlist)
    {
        auto iteR = std::find(Rlist.cbegin(), Rlist.cend(), R);
        auto iR = std::distance(Rlist.cbegin(), iteR);
        for (const auto &Mu_NuqV: coulmat_k)
        {
            const auto Mu = Mu_NuqV.first;
            const int n_mu = basis_abf[Mu];
            for (const auto &Nu_qV: Mu_NuqV.second)
            {
                const auto Nu = Nu_qV.first;
                const int n_nu = basis_abf[Nu];
                coulmat_R[Mu][Nu][R] = std::make_shared<matrix>();
                // a temporary complex matrix to save the transformed matrix
                ComplexMatrix VR_cplx(n_mu, n_nu);
                for (const auto &q_V: Nu_qV.second)
                {
                    auto q = q_V.first;
                    for (auto q_bz: map_irk_ks.at(q))
                    {
                        double ang = - q_bz * (R * latvec) * TWO_PI;
                        complex<double> kphase = complex<double>(cos(ang), sin(ang)) / double(n_k_points);
                        // FIXME: currently support inverse symmetry only
                        if (q_bz == q)
                        {
                            // cout << "Direct:  " << q_bz << " => " << q << ", phase = " << kphase << endl;
                            VR_cplx += (*q_V.second) * kphase;
                        }
                        else
                        {
                            // cout << "Inverse: " << q_bz << " => " << q << ", phase = " << kphase << endl;
                            VR_cplx += conj(*q_V.second) * kphase;
                        }
                    }
                    // minyez debug: check hermicity of Vq
                    // if (iR == 0)
                    // {
                    //     int iq = std::distance(klist.begin(), std::find(klist.begin(), klist.end(), q));
                    //     sprintf(fn, "Vq_Mu_%zu_Nu_%zu_iq_%d.mtx", Mu, Nu, iq);
                    //     print_complex_matrix_mm(*q_V.second, fn);
                    // }
                    // end minyez debug
                }
                *coulmat_R[Mu][Nu][R] = VR_cplx.real();
                // debug print
                // sprintf(fn, "VR_cplx_Mu_%zu_Nu_%zu_iR_%zu.mtx", Mu, Nu, iR);
                // print_complex_matrix_mm(VR_cplx, fn);
                // sprintf(fn, "VR_Mu_%zu_Nu_%zu_iR_%zu.mtx", Mu, Nu, iR);
                // print_matrix_mm(*coulmat_R[Mu][Nu][R], fn);

                // when ordered atom pair is requested, check whether it is available in the original map
                if (return_ordered_atom_pair && Mu != Nu && (coulmat_k.count(Nu) == 0 || coulmat_k.at(Nu).count(Mu) == 0))
                {
                    coulmat_R[Nu][Mu][R] = std::make_shared<matrix>();
                    ComplexMatrix VR_cplx(n_nu, n_mu);
                    for (const auto &q_V: Nu_qV.second)
                    {
                        auto q = q_V.first;
                        for (auto q_bz: map_irk_ks.at(q))
                        {
                            double ang = - q_bz * (R * latvec) * TWO_PI;
                            complex<double> kphase = complex<double>(cos(ang), sin(ang)) / double(n_k_points);
                            if (q_bz == q)
                            {
                                VR_cplx += transpose(*q_V.second, true) * kphase;
                            }
                            else
                            {
                                VR_cplx += transpose(*q_V.second, false) * kphase;
                            }
                        }
                    }
                    *coulmat_R[Nu][Mu][R] = VR_cplx.real();
                }
            }
        }
    }
    // myz debug: check the imaginary part of the coulomb matrix
    // char fn[80];
    /* for (const auto & Mu_NuRV: VR) */
    /* { */
    /*     auto Mu = Mu_NuRV.first; */
    /*     const int n_mu = atom_mu[Mu]; */
    /*     for (const auto & Nu_RV: Mu_NuRV.second) */
    /*     { */
    /*         auto Nu = Nu_RV.first; */
    /*         const int n_nu = atom_mu[Nu]; */
    /*         for (const auto & R_V: Nu_RV.second) */
    /*         { */
    /*             auto R = R_V.first; */
    /*             auto &V = R_V.second; */
    /*             auto iteR = std::find(Rlist.cbegin(), Rlist.cend(), R); */
    /*             auto iR = std::distance(Rlist.cbegin(), iteR); */
    /*             sprintf(fn, "VR_Mu_%zu_Nu_%zu_iR_%zu.mtx", Mu, Nu, iR); */
    /*             print_complex_matrix_mm(*V, fn); */
    /*         } */
    /*     } */
    /* } */
    return coulmat_R;
}

}
