#include "../math/symmetry.h"

#include "testutils.h"
#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <map>
#include <stdexcept>
#include <vector>

using namespace librpa_int;

static Matrix3 lih_primitive_lattice()
{
    return Matrix3(0.0, 0.5, 0.5,
                   0.5, 0.0, 0.5,
                   0.5, 0.5, 0.0);
}

static void add_lih_basis_preserving_fractional_symmetry_operations(
    SpaceGroupSymOps& operations,
    const Matrix3& lattice)
{
    const std::array<std::array<int, 3>, 6> permutations{{
        {{0, 1, 2}},
        {{0, 2, 1}},
        {{1, 0, 2}},
        {{1, 2, 0}},
        {{2, 0, 1}},
        {{2, 1, 0}},
    }};
    const std::array<int, 2> signs{{-1, 1}};

    for (const auto& permutation : permutations)
    {
        for (const int sx : signs)
        {
            for (const int sy : signs)
            {
                for (const int sz : signs)
                {
                    const std::array<int, 3> sign{{sx, sy, sz}};
                    std::array<int, 9> rotation{};
                    for (int col = 0; col != 3; ++col)
                    {
                        rotation[3 * permutation[col] + col] = sign[col];
                    }

                    const Matrix3 col_cartesian_rotation(
                        rotation[0], rotation[1], rotation[2],
                        rotation[3], rotation[4], rotation[5],
                        rotation[6], rotation[7], rotation[8]);

                    SpaceGroupSymOp op;
                    op.rotation =
                        lattice * col_cartesian_rotation.Transpose() * lattice.Inverse();
                    op.translation = {0.0, 0.0, 0.0};
                    op.use_row_convention = true;
                    const Vector3_Order<double> h_position{-0.25, -0.25, -0.25};
                    if (!nearly_integer_vector(
                            multiply_row_vector(h_position, op.rotation) - h_position, 1e-8))
                    {
                        continue;
                    }
                    operations.push_back(op);
                }
            }
        }
    }
    assert(operations.size() == 24);
}

static SpaceGroupSymOps lih_basis_preserving_fractional_symmetry_operations()
{
    SpaceGroupSymOps operations;
    add_lih_basis_preserving_fractional_symmetry_operations(operations, lih_primitive_lattice());
    return operations;
}

static SpaceGroupSymOps add_time_reversal(
    const SpaceGroupSymOps& operations)
{
    SpaceGroupSymOps operations_with_trs;
    operations_with_trs.reserve(operations.size() * 2);
    for (const auto& operation : operations)
    {
        operations_with_trs.push_back(operation);
    }
    for (const auto& operation : operations)
    {
        SpaceGroupSymOp trs_operation = operation;
        trs_operation.rotation = operation.rotation * -1.0;
        operations_with_trs.push_back(trs_operation);
    }
    return operations_with_trs;
}

static double centered_mesh_coordinate(const int i, const int n)
{
    double value = static_cast<double>(i) / static_cast<double>(n);
    if (value > 0.5)
    {
        value -= 1.0;
    }
    return value;
}

static std::vector<Vector3_Order<double>> make_centered_k_grid(const int n)
{
    std::vector<Vector3_Order<double>> kpoints;
    kpoints.reserve(static_cast<std::size_t>(n * n * n));
    for (int ix = 0; ix != n; ++ix)
    {
        for (int iy = 0; iy != n; ++iy)
        {
            for (int iz = 0; iz != n; ++iz)
            {
                kpoints.push_back({centered_mesh_coordinate(ix, n),
                                   centered_mesh_coordinate(iy, n),
                                   centered_mesh_coordinate(iz, n)});
            }
        }
    }
    return kpoints;
}

static std::vector<Vector3_Order<int>> make_centered_R_grid(const int n)
{
    std::vector<Vector3_Order<int>> Rlist;
    Rlist.reserve(static_cast<std::size_t>(n * n * n));
    for (int ix = -n / 2; ix <= (n - 1) / 2; ++ix)
    {
        for (int iy = -n / 2; iy <= (n - 1) / 2; ++iy)
        {
            for (int iz = -n / 2; iz <= (n - 1) / 2; ++iz)
            {
                Rlist.push_back({ix, iy, iz});
            }
        }
    }
    return Rlist;
}

static void assert_kpoint_representatives_match(
    const std::vector<KPointStar>& stars,
    const std::vector<Vector3_Order<double>>& expected_representatives)
{
    assert(stars.size() == expected_representatives.size());
    std::vector<bool> matched(stars.size(), false);
    for (const auto& expected : expected_representatives)
    {
        bool found = false;
        for (std::size_t istar = 0; istar != stars.size(); ++istar)
        {
            const auto& star = stars[istar];
            const auto& representative =
                star.members[static_cast<std::size_t>(star.representative_k_index)].kpoint;
            if (!matched[istar] && same_fractional_kpoint(representative, expected))
            {
                matched[istar] = true;
                found = true;
                break;
            }
        }
        assert(found);
    }
}

static std::size_t count_rspace_sector_entries(const SpaceGroupRSpaceSector& sector)
{
    std::size_t count = 0;
    for (const auto& pair_Rs : sector)
    {
        count += pair_Rs.second.size();
    }
    return count;
}

static void test_rotation()
{
    Vector3_Order<double> a1{0.3, 0.2, 0.1};
    const Matrix3 rotation(0.0, -1.0, 0.0,
                           1.0, 0.0, 0.0,
                           0.0, 0.0, 1.0);
    SpaceGroupSymOp op1{rotation, {0.0, 0.0, 0.0}};
    auto a2 = apply_space_group_symmetry_operation(op1, a1);
    assert(fequal(a2.x, 0.2));
    assert(fequal(a2.y, -0.3));
    assert(fequal(a2.z, 0.1));
}

static void test_composition()
{
    const Matrix3 lhs_rotation(0.0, -1.0, 0.0,
                               1.0, 0.0, 0.0,
                               0.0, 0.0, 1.0);
    const Matrix3 rhs_rotation(1.0, 0.0, 0.0,
                               0.0, 0.0, -1.0,
                               0.0, 1.0, 0.0);
    const SpaceGroupSymOp lhs{lhs_rotation, {0.1, 0.2, 0.3}};
    const SpaceGroupSymOp rhs{rhs_rotation, {0.4, 0.5, 0.6}};
    const auto composed = compose_space_group_symmetry_operations(lhs, rhs);

    const Vector3_Order<double> coord{0.7, 0.8, 0.9};
    const auto sequential =
        apply_space_group_symmetry_operation(rhs,
                                             apply_space_group_symmetry_operation(lhs, coord));
    const auto direct = apply_space_group_symmetry_operation(composed, coord);

    assert(fequal(direct.x, sequential.x));
    assert(fequal(direct.y, sequential.y));
    assert(fequal(direct.z, sequential.z));
}

static void test_column_convention()
{
    const Matrix3 lhs_rotation(0.0, -1.0, 0.0,
                               1.0, 0.0, 0.0,
                               0.0, 0.0, 1.0);
    const Matrix3 rhs_rotation(1.0, 0.0, 0.0,
                               0.0, 0.0, -1.0,
                               0.0, 1.0, 0.0);
    const SpaceGroupSymOp lhs{lhs_rotation, {0.1, 0.2, 0.3}, false};
    const SpaceGroupSymOp rhs{rhs_rotation, {0.4, 0.5, 0.6}, false};

    const Vector3_Order<double> coord{0.7, 0.8, 0.9};
    const auto applied = apply_space_group_symmetry_operation(lhs, coord);
    assert(fequal(applied.x, -0.7));
    assert(fequal(applied.y, 0.9));
    assert(fequal(applied.z, 1.2));

    const auto composed = compose_space_group_symmetry_operations(lhs, rhs);
    const auto sequential =
        apply_space_group_symmetry_operation(rhs,
                                             apply_space_group_symmetry_operation(lhs, coord));
    const auto direct = apply_space_group_symmetry_operation(composed, coord);
    assert(fequal(direct.x, sequential.x));
    assert(fequal(direct.y, sequential.y));
    assert(fequal(direct.z, sequential.z));
}

static void test_row_fractional_rotation_to_cartesian()
{
    const Matrix3 lattice_vectors(2.0, 0.0, 0.0,
                                  0.5, 1.5, 0.0,
                                  0.2, 0.3, 2.0);
    const Matrix3 rotation(0.0, -1.0, 0.0,
                           1.0, 0.0, 0.0,
                           0.0, 0.0, 1.0);
    const Vector3_Order<double> frac{0.2, 0.3, 0.4};

    const auto cartesian_rotation =
        row_fractional_rotation_to_cartesian(rotation, lattice_vectors);
    const auto cart = Vector3_Order<double>(lattice_vectors.Transpose() * frac);
    const auto rotated_frac = multiply_row_vector(frac, rotation);
    const auto rotated_cart = Vector3_Order<double>(lattice_vectors.Transpose() * rotated_frac);
    const auto direct_cart = Vector3_Order<double>(cartesian_rotation * cart);

    assert(fequal(direct_cart.x, rotated_cart.x));
    assert(fequal(direct_cart.y, rotated_cart.y));
    assert(fequal(direct_cart.z, rotated_cart.z));
}

static void test_rounded_hexagonal_rspace_sector_is_scale_invariant()
{
    const Matrix3 c3(-1.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0);
    const SpaceGroupSymOps operations{SpaceGroupSymOp{Matrix3::IDENTITY, {0.0, 0.0, 0.0}, true},
                                      SpaceGroupSymOp{c3, {0.0, 0.0, 0.0}, true},
                                      SpaceGroupSymOp{c3 * c3, {0.0, 0.0, 0.0}, true}};
    const std::map<int, Vector3_Order<double>> positions{{0, {0.0, 0.0, 0.0}}};
    const std::map<int, int> types{{0, 0}};
    // The origin and the six hexagonal neighbours form three C3 orbits.
    const std::vector<Vector3_Order<int>> neighbours{{0, 0, 0},  {1, 0, 0},  {0, 1, 0}, {-1, 1, 0},
                                                     {-1, 0, 0}, {0, -1, 0}, {1, -1, 0}};
    for (const double scale : {0.1, 1.0, 10.0})
    {
        // sqrt(3)/2 is rounded at realistic input precision. An absolute
        // lattice-metric threshold wrongly loses C3 when the cell is scaled.
        const Matrix3 lattice = Matrix3(10.0, 0.0, 0.0, 5.0, 8.660254, 0.0, 0.0, 0.0, 13.0) * scale;
        const auto sector = build_space_group_rspace_irreducible_sector(
            operations, positions, types, neighbours, &lattice);
        assert(count_rspace_sector_entries(sector) == 3);
    }
}

static void test_atom_to_inequivalent_symmetry_mapping()
{
    const Matrix3 lattice_vectors;
    const std::vector<Vector3_Order<double>> atom_positions{
        {0.0, 0.0, 0.0},
        {0.5, 0.5, 0.0},
        {0.25, 0.25, 0.0},
        {0.75, 0.75, 0.0},
    };

    SpaceGroupSymOps operations;
    operations.push_back(SpaceGroupSymOp{Matrix3(), {0.0, 0.0, 0.0}});
    operations.push_back(SpaceGroupSymOp{Matrix3(), {0.5, 0.5, 0.0}});

    const auto mappings = build_atom_to_inequivalent_symmetry_mapping(
        atom_positions, lattice_vectors, operations);
    const auto inequivalent_atoms = collect_inequivalent_atoms(mappings);

    assert(inequivalent_atoms.size() == 2);
    assert(inequivalent_atoms[0] == 0);
    assert(inequivalent_atoms[1] == 2);

    assert(mappings[0].inequivalent_atom == 0);
    assert(mappings[0].isym == 0);
    assert(mappings[0].return_lattice == Vector3_Order<int>(0, 0, 0));

    assert(mappings[1].inequivalent_atom == 0);
    assert(mappings[1].isym == 1);
    assert(mappings[1].return_lattice == Vector3_Order<int>(1, 1, 0));

    assert(mappings[2].inequivalent_atom == 2);
    assert(mappings[2].isym == 0);

    assert(mappings[3].inequivalent_atom == 2);
    assert(mappings[3].isym == 1);
    assert(mappings[3].return_lattice == Vector3_Order<int>(1, 1, 0));
}

static void test_symmetry_mapping_keeps_operations_fractional()
{
    const Matrix3 lattice_vectors(2.0, 0.0, 0.0, 0.5, 1.5, 0.0, 0.2, 0.3, 2.0);
    const Matrix3 rotation(0.0, -1.0, 0.0,
                           1.0, 0.0, 0.0,
                           0.0, 0.0, 1.0);
    const Vector3_Order<double> source_frac{0.2, 0.3, 0.4};
    const Vector3_Order<double> translation_frac{0.25, 0.1, 0.0};
    const Vector3_Order<double> target = restrict_fractional_coordinate(
        multiply_row_vector(source_frac, rotation) + translation_frac, 1e-10);

    const std::vector<Vector3_Order<double>> atom_positions{target, source_frac};
    SpaceGroupSymOps fractional_operations;
    fractional_operations.push_back(SpaceGroupSymOp{Matrix3(), {0.0, 0.0, 0.0}});
    fractional_operations.push_back(SpaceGroupSymOp{rotation, translation_frac});

    const auto mappings = build_atom_to_inequivalent_symmetry_mapping(
        atom_positions, lattice_vectors, fractional_operations, 1e-10);
    const auto inequivalent_atoms = collect_inequivalent_atoms(mappings);

    assert(inequivalent_atoms.size() == 1);
    assert(inequivalent_atoms[0] == 0);
    assert(mappings[1].inequivalent_atom == 0);
    assert(mappings[1].isym == 1);
    assert(mappings[1].return_lattice == Vector3_Order<int>(0, -1, 0));
}

static void test_kpoint_rotation_and_target_fold()
{
    const Matrix3 rotation(1.0, 1.0, 0.0,
                           0.0, 1.0, 0.0,
                           0.0, 0.0, 1.0);
    const SpaceGroupSymOp op{rotation, {0.0, 0.0, 0.0}};
    const Vector3_Order<double> rpoint{0.2, 0.3, 0.0};
    const Vector3_Order<double> kpoint{0.4, 0.6, 0.0};

    const auto rotated_r = apply_space_group_symmetry_operation(op, rpoint);
    const auto rotated_k = apply_space_group_rotation_to_kpoint(op, kpoint);
    assert(fequal(kpoint * rpoint, rotated_k * rotated_r));
    assert(fequal(rotated_k.x, -0.2));
    assert(fequal(rotated_k.y, 0.6));
    assert(fequal(rotated_k.z, 0.0));

    const std::vector<Vector3_Order<double>> target_kpoints{
        {1.8, 0.6, 0.0},
        {0.4, 0.6, 0.0},
    };
    const auto target_folded = fold_fractional_kpoint_to_targets(rotated_k, target_kpoints);
    assert(target_folded.target_k_index == 0);
    assert(target_folded.kpoint == Vector3_Order<double>(1.8, 0.6, 0.0));
    assert(target_folded.fold_G == Vector3_Order<int>(-2, 0, 0));

    assert(same_fractional_kpoint({0.0, 0.5, -1.0}, {1.0, -0.5, 0.0}));
    FoldedKPoint folded_to_single_target;
    assert(try_fold_fractional_kpoint_to_target(
        {1.25, -0.75, 0.0}, {0.25, 0.25, 0.0}, 1e-8, folded_to_single_target));
    assert(folded_to_single_target.kpoint == Vector3_Order<double>(0.25, 0.25, 0.0));
    assert(folded_to_single_target.fold_G == Vector3_Order<int>(1, -1, 0));
}

static void test_kpoint_stars_from_full_grid()
{
    const Matrix3 identity;
    // Clockwise C4 rotations; check x' = x R
    const Matrix3 rot90(0.0, -1.0, 0.0,
                        1.0, 0.0, 0.0,
                        0.0, 0.0, 1.0);
    const Matrix3 rot180(-1.0, 0.0, 0.0,
                         0.0, -1.0, 0.0,
                         0.0, 0.0, 1.0);
    const Matrix3 rot270(0.0, 1.0, 0.0,
                         -1.0, 0.0, 0.0,
                         0.0, 0.0, 1.0);
    SpaceGroupSymOps operations;
    operations.push_back(SpaceGroupSymOp{identity, {0.0, 0.0, 0.0}});
    operations.push_back(SpaceGroupSymOp{rot90, {0.0, 0.0, 0.0}});
    operations.push_back(SpaceGroupSymOp{rot180, {0.0, 0.0, 0.0}});
    operations.push_back(SpaceGroupSymOp{rot270, {0.0, 0.0, 0.0}});

    const std::vector<Vector3_Order<double>> full_kpoints{
        {0.0, 0.0, 0.0},
        {0.5, 0.0, 0.0},
        {0.0, -0.5, 0.0},
        {0.5, 0.5, 0.0},
    };

    const auto rotate_k = apply_space_group_rotation_to_kpoint(operations[1], full_kpoints[1]);
    assert(fequal(rotate_k.x, 0.0));
    assert(fequal(rotate_k.y, -0.5));
    assert(fequal(rotate_k.z, 0.0));

    const auto stars = build_kpoint_stars(full_kpoints, operations);
    assert(stars.size() == 3);
    assert(stars[0].members.size() == 1);
    assert(stars[0].sym_mappings.size() == stars[0].members.size());
    assert(stars[0].representative_k_index == 0);
    assert(stars[0].members[stars[0].representative_k_index].full_k_index == 0);
    assert(stars[1].members.size() == 2);
    assert(stars[1].sym_mappings.size() == stars[1].members.size());
    assert(stars[1].representative_k_index == 0);
    assert(stars[1].members[stars[1].representative_k_index].full_k_index == 1);
    assert(stars[1].members[0].full_k_index == 1);
    assert(stars[1].sym_mappings[0].isym == 0);
    assert(stars[1].sym_mappings[0].fold_G == Vector3_Order<int>(0, 0, 0));
    assert(stars[1].members[1].full_k_index == 2);
    assert(stars[1].members[1].kpoint == Vector3_Order<double>(0.0, -0.5, 0.0));
    assert(stars[1].sym_mappings[1].isym == 1);
    assert(stars[1].sym_mappings[1].fold_G == Vector3_Order<int>(0, 0, 0));
    assert(stars[2].members.size() == 1);
    assert(stars[2].sym_mappings.size() == stars[2].members.size());
    assert(stars[2].representative_k_index == 0);
    assert(stars[2].members[stars[2].representative_k_index].full_k_index == 3);

    const std::vector<Vector3_Order<double>> preferred_representatives{{0.0, -0.5, 0.0}};
    const auto hinted_stars = build_kpoint_stars(full_kpoints, operations, preferred_representatives);
    assert(hinted_stars[1].members.size() == 2);
    assert(hinted_stars[1].representative_k_index == 1);
    assert(hinted_stars[1].members[hinted_stars[1].representative_k_index].full_k_index == 2);
    assert(hinted_stars[1].sym_mappings[0].isym == 1);
    assert(hinted_stars[1].sym_mappings[0].fold_G == Vector3_Order<int>(-1, 0, 0));
    assert(hinted_stars[1].sym_mappings[1].isym == 0);
    assert(hinted_stars[1].sym_mappings[1].fold_G == Vector3_Order<int>(0, 0, 0));
}

static void test_lih_rocksalt_kstars_match_abacus()
{
    const auto operations =
        add_time_reversal(lih_basis_preserving_fractional_symmetry_operations());

    {
        const std::vector<Vector3_Order<double>> expected_representatives{
            {0.0, 0.0, 0.0},
            {0.25, 0.25, 0.25},
            {0.5, 0.5, 0.5},
            {0.25, 0.25, 0.0},
            {0.5, 0.5, 0.25},
            {0.5, 0.25, 0.25},
            {0.5, 0.5, 0.0},
            {-0.25, 0.5, 0.25},
        };
        const auto stars =
            build_kpoint_stars(make_centered_k_grid(4), operations, expected_representatives);
        assert_kpoint_representatives_match(stars, expected_representatives);
        std::size_t member_count = 0;
        for (const auto& star : stars)
        {
            member_count += star.members.size();
        }
        assert(member_count == 64);
    }

    {
        const std::vector<Vector3_Order<double>> expected_representatives{
            {0.0, 0.0, 0.0},
            {0.2, 0.2, 0.2},
            {0.4, 0.4, 0.4},
            {0.2, 0.2, 0.0},
            {0.4, 0.4, 0.2},
            {-0.4, 0.4, 0.4},
            {0.2, 0.4, 0.2},
            {0.4, 0.4, 0.0},
            {-0.2, 0.4, 0.4},
            {-0.4, 0.4, 0.2},
        };
        const auto stars =
            build_kpoint_stars(make_centered_k_grid(5), operations, expected_representatives);
        assert_kpoint_representatives_match(stars, expected_representatives);
        std::size_t member_count = 0;
        for (const auto& star : stars)
        {
            member_count += star.members.size();
        }
        assert(member_count == 125);
    }

    {
        const std::vector<Vector3_Order<double>> expected_representatives{
            {0.0, 0.0, 0.0},
            {1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0},
            {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0},
            {0.5, 0.5, 0.5},
            {1.0 / 6.0, 1.0 / 6.0, 0.0},
            {1.0 / 3.0, 1.0 / 6.0, 1.0 / 3.0},
            {0.5, 0.5, 1.0 / 3.0},
            {0.5, 1.0 / 3.0, 1.0 / 3.0},
            {1.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0},
            {1.0 / 3.0, 1.0 / 3.0, 0.0},
            {0.5, 0.5, 1.0 / 6.0},
            {-1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0},
            {0.5, 0.5, 0.0},
            {0.5, 1.0 / 3.0, 1.0 / 6.0},
            {-1.0 / 3.0, 0.5, 1.0 / 3.0},
            {-1.0 / 6.0, 0.5, 1.0 / 3.0},
        };
        const auto stars =
            build_kpoint_stars(make_centered_k_grid(6), operations, expected_representatives);
        assert_kpoint_representatives_match(stars, expected_representatives);
        std::size_t member_count = 0;
        for (const auto& star : stars)
        {
            member_count += star.members.size();
        }
        assert(member_count == 216);
    }
}

static void assert_lih_rspace_sector_uses_input_fractional_representatives(
    const int n,
    const std::array<std::size_t, 4>& expected_pair_counts,
    const std::size_t expected_total_count)
{
    const auto operations = lih_basis_preserving_fractional_symmetry_operations();
    const std::map<int, Vector3_Order<double>> coord_frac{
        {0, {0.0, 0.0, 0.0}},
        {1, {-0.25, -0.25, -0.25}},
    };
    const std::map<int, int> atom_to_type{{0, 0}, {1, 1}};
    const auto lattice = lih_primitive_lattice();

    const auto sector = build_space_group_rspace_irreducible_sector(
        operations, coord_frac, atom_to_type, make_centered_R_grid(n), &lattice);

    assert(sector.at({0, 0}).size() == expected_pair_counts[0]);
    assert(sector.at({0, 1}).size() == expected_pair_counts[1]);
    assert(sector.at({1, 0}).size() == expected_pair_counts[2]);
    assert(sector.at({1, 1}).size() == expected_pair_counts[3]);
    assert(count_rspace_sector_entries(sector) == expected_total_count);
}

static void test_lih_rocksalt_rspace_irreducible_sectors_use_input_fractional_representatives()
{
    assert_lih_rspace_sector_uses_input_fractional_representatives(4, {{15, 16, 13, 15}}, 59);
    assert_lih_rspace_sector_uses_input_fractional_representatives(5, {{22, 22, 22, 22}}, 88);
    assert_lih_rspace_sector_uses_input_fractional_representatives(6, {{37, 39, 34, 37}}, 147);
}

static void test_get_space_group_atom_mapping_mgo()
{
    // MgO FCC
    const std::map<size_t, Vector3_Order<double>> coord_frac{{0, {0.0, 0.0, 0.0}},
                                                             {1, {0.5, 0.5, 0.5}}};
    const std::map<size_t, int> atom_to_type{{0, 0}, {1, 1}};
    const auto map_iden = get_space_group_atom_mapping(SpaceGroupSymOp::IDENTITY, coord_frac, atom_to_type);
    assert(map_iden.atom_map[0] == 0);
    assert(map_iden.atom_map[1] == 1);
    {
        const auto map_inv = get_space_group_atom_mapping(SpaceGroupSymOp::INVERSE, coord_frac, atom_to_type);
        assert(map_inv.atom_map[0] == 0);
        assert(map_inv.atom_map[1] == 1);
        const auto return_lattice = Vector3_Order<int>{-1, -1, -1};
        std::cout << "map_inv.return_lattice[1] " << map_inv.return_lattice[1] << std::endl;
        assert(map_inv.return_lattice[1] == return_lattice);
    }
    {
        const auto map_c41z = get_space_group_atom_mapping(SpaceGroupSymOp::C41_Z, coord_frac, atom_to_type);
        assert(map_c41z.atom_map[0] == 0);
        assert(map_c41z.atom_map[1] == 1);
        const auto return_lattice = Vector3_Order<int>{-1, 0, 0};
        std::cout << "map_c41z.return_lattice[1] " << map_c41z.return_lattice[1] << std::endl;
        assert(map_c41z.return_lattice[1] == return_lattice);
    }
    {
        SpaceGroupSymOp op;
        op.rotation = {1, 0, 0, 0, 0, -1, 0, 1, 0};
        op.use_row_convention = true;
        const auto map_row = get_space_group_atom_mapping(op, coord_frac, atom_to_type);
        auto return_lattice = Vector3_Order<int>{0, 0, -1};
        std::cout << "map_row.return_lattice[1] " << map_row.return_lattice[1] << std::endl;
        assert(map_row.atom_map[0] == 0);
        assert(map_row.atom_map[1] == 1);
        assert(map_row.return_lattice[1] == return_lattice);
        op.use_row_convention = false;
        const auto map_col = get_space_group_atom_mapping(op, coord_frac, atom_to_type);
        std::cout << "map_col.return_lattice[1] " << map_col.return_lattice[1] << std::endl;
        return_lattice = Vector3_Order<int>{0, -1, 0};
        assert(map_col.atom_map[0] == 0);
        assert(map_col.atom_map[1] == 1);
        assert(map_col.return_lattice[1] == return_lattice);
    }
}

static void test_return_lattice_preserves_input_fractional_representative()
{
    const std::map<size_t, Vector3_Order<double>> coord_frac{{0, {0.0, 0.0, 0.0}},
                                                             {1, {-0.25, -0.25, -0.25}}};
    const std::map<size_t, int> atom_to_type{{0, 0}, {1, 1}};
    SpaceGroupSymOp op;
    op.rotation = {0, 0, -1, 1, 0, -1, 0, 1, -1};
    op.use_row_convention = true;

    const auto map = get_space_group_atom_mapping(op, coord_frac, atom_to_type);
    assert(map.atom_map[1] == 1);
    assert(map.return_lattice[1] == Vector3_Order<int>(0, 0, 1));
}

static void test_direct_fractional_coordinates_fix_si_symmetry_atom_mapping()
{
    const std::map<size_t, int> atom_to_type{{0, 0}, {1, 0}};
    SpaceGroupSymOp op;
    op.rotation = {0, -1, 0, 1, -1, 0, 0, -1, 1};
    op.translation = {0.0, 0.5, -1.0};
    op.use_row_convention = true;

    const std::map<size_t, Vector3_Order<double>> coord_frac_from_mismatched_cart{
        {0, {0.12500258782275, 0.12500258782275, 0.12500258782275}},
        {1, {0.87501811475925, 0.87501811475925, 0.87501811475925}},
    };
    bool failed_with_mismatched_cart = false;
    try
    {
        (void)get_space_group_atom_mapping(
            op, coord_frac_from_mismatched_cart, atom_to_type, 5e-5);
    }
    catch (const std::runtime_error&)
    {
        failed_with_mismatched_cart = true;
    }
    assert(failed_with_mismatched_cart);

    const std::map<size_t, Vector3_Order<double>> coord_frac_from_direct_stru{
        {0, {0.125, 0.125, 0.125}},
        {1, {0.875, 0.875, 0.875}},
    };
    const auto map =
        get_space_group_atom_mapping(op, coord_frac_from_direct_stru, atom_to_type, 5e-5);
    assert(map.atom_map[0] == 0);
    assert(map.atom_map[1] == 1);
}

int main()
{
    test_rotation();
    test_composition();
    test_column_convention();
    test_row_fractional_rotation_to_cartesian();
    test_atom_to_inequivalent_symmetry_mapping();
    test_symmetry_mapping_keeps_operations_fractional();
    test_kpoint_rotation_and_target_fold();
    test_kpoint_stars_from_full_grid();
    test_lih_rocksalt_kstars_match_abacus();
    test_lih_rocksalt_rspace_irreducible_sectors_use_input_fractional_representatives();
    test_rounded_hexagonal_rspace_sector_is_scale_invariant();
    test_get_space_group_atom_mapping_mgo();
    test_return_lattice_preserves_input_fractional_representative();
    test_direct_fractional_coordinates_fix_si_symmetry_atom_mapping();
    return 0;
}
