#include <cmath>
#include <fstream>
#include <iostream>

#include "../core/pbc.h"
#include "../core/symmetry_context.h"
using namespace librpa_int;
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
        double err = 0., dense_err = 0., wfc_err = 0., norm = 0.;
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
    }
    return failed ? 1 : 0;
}
