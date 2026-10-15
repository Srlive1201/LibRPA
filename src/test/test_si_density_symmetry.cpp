#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>

#include "../core/pbc.h"
#include "../core/symmetry_context.h"
#include "../utils/constants.h"
using namespace librpa_int;

static SymmetryKStarMember generate_reference_member(
    SymmetryContext& ctx, PeriodicBoundaryData pbc,
    const Vector3_Order<double>& k_ibz, const Vector3_Order<double>& k_bz)
{
    // Close the fixture operation into a cyclic subgroup. The reference density
    // remains producer data; only the route to that density is generated here.
    const auto generator = ctx.rspace_operations.at(1);
    ctx.rspace_operations = {SpaceGroupSymOp::IDENTITY};
    auto operation = generator;
    for (int order = 1; !operation.is_identity(); ++order)
    {
        if (order > 48) throw std::runtime_error("Fixture subgroup does not close");
        ctx.rspace_operations.push_back(operation);
        operation = compose_space_group_symmetry_operations(operation, generator);
        operation.translation = restrict_fractional_coordinate(operation.translation);
    }
    auto operations_with_tr = ctx.rspace_operations;
    for (auto op : ctx.rspace_operations)
    {
        op.rotation = op.rotation * -1.0;
        operations_with_tr.push_back(op);
    }
    const auto full_kpoints = build_uniform_kmesh_frac({4, 4, 4});
    std::vector<double> kvecs;
    for (const auto& kfrac : full_kpoints)
    {
        const auto kcart = kfrac * pbc.G;
        kvecs.insert(kvecs.end(), {kcart.x * TWO_PI, kcart.y * TWO_PI, kcart.z * TWO_PI});
    }
    pbc.set_kgrids_kvec(4, 4, 4, kvecs);
    const auto stars = build_kpoint_stars(full_kpoints, operations_with_tr, {k_ibz}, 1e-5);
    pbc.kfrac_list.clear();
    pbc.klist.clear();
    pbc.klist_coul.clear();
    for (const auto& star : stars)
    {
        const auto& representative = star.members.at(star.representative_k_index).kpoint;
        pbc.kfrac_list.push_back(representative);
        pbc.klist.push_back(representative * pbc.G);
        pbc.klist_coul.push_back(representative * pbc.G);
    }
    pbc.map_irk_ks.clear();
    ctx.build_periodic_mappings(pbc, {});
    ctx.build_kstar_member_rotations(2);
    for (const auto& star : ctx.kstars)
    {
        if (!same_fractional_kpoint(star.k_ibz, k_ibz)) continue;
        for (const auto& member : star.members)
            if (same_fractional_kpoint(member.k_bz, k_bz)) return member;
    }
    throw std::runtime_error("Automatically generated star lost the reference member");
}

int main()
{
    std::ifstream in(SI_DENSITY_REFERENCE);
    if (!in)
    {
        std::cerr << "Cannot open Si density reference" << std::endl;
        return 2;
    }
    int n = 0;
    if (!(in >> n) || n != 3)
    {
        std::cerr << "Expected three Si density reference cases" << std::endl;
        return 2;
    }
    bool failed = false;

    // Convention guard, independent of the density reference file below.
    // In the row-major convention C_bz = C_ibz * rotation the block of the
    // source atom I must be placed at the rows of S(I) and the columns of I.
    // This is checked through the return code rather than assert(), so a
    // regression of the placement (or a transposed placement) fails the test
    // in every build type, including Release with NDEBUG.
    {
        SymmetryContext ctx;
        ctx.rspace_operations.push_back(SpaceGroupSymOp::IDENTITY);
        SymmetryKStarMember member;
        member.spatial_isym = 0;
        member.k_bz = {0.0, 0.0, 0.0};
        SymmetryKAtomRotation rot_0;
        rot_0.atom_from = 0;
        rot_0.atom_to = 1;
        rot_0.atom_type = 0;
        rot_0.lmax = 0;
        rot_0.bloch_rsh_rotations[0] = ComplexMatrix(1, 1);
        rot_0.bloch_rsh_rotations[0](0, 0) = {2.0, 0.5};
        SymmetryKAtomRotation rot_1;
        rot_1.atom_from = 1;
        rot_1.atom_to = 0;
        rot_1.atom_type = 0;
        rot_1.lmax = 0;
        rot_1.bloch_rsh_rotations[0] = ComplexMatrix(1, 1);
        rot_1.bloch_rsh_rotations[0](0, 0) = {3.0, -0.25};
        member.atom_rotations = {rot_0, rot_1};
        SpeciesBasisLayout layout;
        layout.label = "X";
        layout.set({0});
        const std::map<atom_t, size_t> atom_nabf{{0, 1}, {1, 1}};
        const auto rotation = build_symmetry_kspace_rotation_matrix(
            ctx, {layout}, member, atom_nabf, {0.0, 0.0, 0.0});
        const std::complex<double> block_0{2.0, 0.5};
        const std::complex<double> block_1{3.0, -0.25};
        const bool placement_ok = std::abs(rotation(1, 0) - block_0) < 1e-12
                                  && std::abs(rotation(0, 1) - block_1) < 1e-12
                                  && std::abs(rotation(0, 0)) < 1e-12
                                  && std::abs(rotation(1, 1)) < 1e-12;
        std::cout << "rotation block placement (rows atom_to, columns atom_from): "
                  << (placement_ok ? "ok" : "FAILED") << std::endl;
        failed |= !placement_ok;
    }

    for (int c = 0; c < n; ++c)
    {
        double r[9] = {}, t[3] = {}, ki[3] = {}, kb[3] = {};
        for (auto& x : r) in >> x;
        for (auto& x : t) in >> x;
        for (auto& x : ki) in >> x;
        for (auto& x : kb) in >> x;
        if (!in)
        {
            std::cerr << "Invalid Si symmetry reference metadata" << std::endl;
            return 2;
        }
        SymmetryContext ctx;
        PeriodicBoundaryData pbc;
        pbc.set_latvec({0., 1., 1., 1., 0., 1., 1., 1., 0.});
        ctx.set_crystal_structure(pbc.latvec, pbc.G, {{0, 0}, {1, 0}},
                                  {{0, {0., 0., 0.}}, {1, {.25, .25, .25}}});
        ctx.basis_convention = {-1, 0, LIBRPA_ANGULAR_ORDER_ABS_PM, LIBRPA_RSH_COEFF_M_1,
                                LIBRPA_RSH_COEFF_1_M};
        ctx.rspace_operations.push_back(SpaceGroupSymOp::IDENTITY);
        ctx.rspace_operations.push_back(SpaceGroupSymOp{
            Matrix3(r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7], r[8]), {t[0], t[1], t[2]}});
        Vector3_Order<double> k_ibz{ki[0], ki[1], ki[2]}, k_bz{kb[0], kb[1], kb[2]};
        auto member = build_symmetry_kspace_operation_member(ctx, 1, false, k_bz, k_ibz, 2);
        SpeciesBasisLayout layout;
        layout.label = "Si";
        layout.set({0, 0, 0, 1, 1, 1, 2, 2});
        symmetry_atom_block_matrix_map_t blocks;
        ComplexMatrix density_ibz(44, 44);
        for (int a = 0; a < 2; ++a)
            for (int b = 0; b < 2; ++b) blocks[a][b] = ComplexMatrix(22, 22);
        for (int i = 0; i < 44; ++i)
            for (int j = 0; j < 44; ++j)
            {
                double re = 0., im = 0.;
                if (!(in >> re >> im))
                {
                    std::cerr << "Invalid IBZ density reference" << std::endl;
                    return 2;
                }
                blocks[i / 22][j / 22](i % 22, j % 22) = {re, im};
                density_ibz(i, j) = {re, im};
            }
        auto actual = rotate_symmetry_kspace_operator_blocks(
            ctx, {layout}, member, blocks, {{0, 22}, {1, 22}}, k_ibz, false, nullptr, &k_bz);
        const auto dense_actual = rotate_symmetry_kspace_matrix(
            ctx, {layout}, member, density_ibz, {{0, 22}, {1, 22}}, k_ibz, false, &k_bz);
        // C_BZ = C_IBZ * rotation, hence D_BZ = rotation^T * D_IBZ * rotation^*.
        // Compare the row-major wave-function route to the independent density too.
        const auto rotation = build_symmetry_kspace_rotation_matrix(
            ctx, {layout}, member, {{0, 22}, {1, 22}}, k_ibz, false, &k_bz);
        const auto wfc_density = transpose(rotation, false) * density_ibz * conj(rotation);
        auto automatic_ctx = ctx;
        const auto automatic_member = generate_reference_member(automatic_ctx, pbc, k_ibz, k_bz);
        const auto automatic_actual = rotate_symmetry_kspace_operator_blocks(
            automatic_ctx, {layout}, automatic_member, blocks, {{0, 22}, {1, 22}},
            k_ibz, false, nullptr, &k_bz);
        const auto automatic_dense = rotate_symmetry_kspace_matrix(
            automatic_ctx, {layout}, automatic_member, density_ibz, {{0, 22}, {1, 22}},
            k_ibz, false, &k_bz);
        const auto automatic_rotation = build_symmetry_kspace_rotation_matrix(
            automatic_ctx, {layout}, automatic_member, {{0, 22}, {1, 22}}, k_ibz, false, &k_bz);
        const auto automatic_wfc = transpose(automatic_rotation, false)
            * (automatic_member.time_reversal ? conj(density_ibz) : density_ibz)
            * conj(automatic_rotation);
        double err = 0., dense_err = 0., wfc_err = 0., norm = 0.;
        double auto_err = 0., auto_dense_err = 0., auto_wfc_err = 0.;
        for (int i = 0; i < 44; ++i)
            for (int j = 0; j < 44; ++j)
            {
                double re = 0., im = 0.;
                if (!(in >> re >> im))
                {
                    std::cerr << "Invalid full-grid density reference" << std::endl;
                    return 2;
                }
                std::complex<double> ref{re, im};
                err += std::norm(actual.at(i / 22).at(j / 22)(i % 22, j % 22) - ref);
                dense_err += std::norm(dense_actual(i, j) - ref);
                wfc_err += std::norm(wfc_density(i, j) - ref);
                auto_err += std::norm(automatic_actual.at(i / 22).at(j / 22)(i % 22, j % 22) - ref);
                auto_dense_err += std::norm(automatic_dense(i, j) - ref);
                auto_wfc_err += std::norm(automatic_wfc(i, j) - ref);
                norm += std::norm(ref);
            }
        if (!in || norm <= 0.)
        {
            std::cerr << "Invalid density reference" << std::endl;
            return 2;
        }
        const double relative = std::sqrt(err / norm);
        const double dense_relative = std::sqrt(dense_err / norm);
        const double wfc_relative = std::sqrt(wfc_err / norm);
        std::cout << "case " << c << " block_relative_error " << relative
                  << " dense_relative_error " << dense_relative << " wfc_relative_error "
                  << wfc_relative << std::endl;
        failed |= !std::isfinite(relative) || relative > 1e-3 || !std::isfinite(dense_relative) ||
                  dense_relative > 1e-3 || !std::isfinite(wfc_relative) || wfc_relative > 1e-3;
        const double auto_relative = std::sqrt(auto_err / norm);
        const double auto_dense_relative = std::sqrt(auto_dense_err / norm);
        const double auto_wfc_relative = std::sqrt(auto_wfc_err / norm);
        std::cout << "case " << c << " automatic_block_relative_error " << auto_relative
                  << " automatic_dense_relative_error " << auto_dense_relative
                  << " automatic_wfc_relative_error " << auto_wfc_relative
                  << " time_reversal " << automatic_member.time_reversal << std::endl;
        failed |= !std::isfinite(auto_relative) || auto_relative > 1e-3
                  || !std::isfinite(auto_dense_relative) || auto_dense_relative > 1e-3
                  || !std::isfinite(auto_wfc_relative) || auto_wfc_relative > 1e-3;
    }
    return failed ? 1 : 0;
}
