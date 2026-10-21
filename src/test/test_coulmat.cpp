#include <cassert>
#include <cmath>
#include <complex>
#include <iostream>
#include <map>
#include <memory>
#include <vector>

#include "../core/coulmat.h"
#include "../io/global_io.h"
#include "../mpi/base_mpi.h"
#include "../utils/constants.h"

using namespace librpa_int;

namespace
{

void add_scalar_vq_block(atpair_k_cplx_mat_t& vq, const atom_t atom_i, const atom_t atom_j,
                         const Vector3_Order<double>& q, const std::complex<double> value)
{
    auto block = std::make_shared<ComplexMatrix>(1, 1);
    (*block)(0, 0) = value;
    vq[atom_i][atom_j][q] = std::move(block);
}

void add_scalar_vq_blocks(atpair_k_cplx_mat_t& vq, const Vector3_Order<double>& q,
                          const std::map<atom_t, std::map<atom_t, std::complex<double>>>& blocks)
{
    for (const auto& [atom_i, row] : blocks)
    {
        for (const auto& [atom_j, value] : row)
        {
            add_scalar_vq_block(vq, atom_i, atom_j, q, value);
        }
    }
}

void scale_replicated_collective_input(atpair_k_cplx_mat_t& vq)
{
    const double scale = 1.0 / static_cast<double>(global::mpi_comm_global_h.nprocs);
    for (auto& [atom_i, row] : vq)
    {
        (void)atom_i;
        for (auto& [atom_j, q_blocks] : row)
        {
            (void)atom_j;
            for (auto& [q, block] : q_blocks)
            {
                (void)q;
                *block *= scale;
            }
        }
    }
}

PeriodicBoundaryData make_full_pbc()
{
    PeriodicBoundaryData pbc;
    pbc.set_latvec({1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0});
    const std::vector<double> kvecs{0.0, 0.0, 0.0, TWO_PI / 3.0, 0.0, 0.0, TWO_PI * 2.0 / 3.0,
                                    0.0, 0.0};
    pbc.set_kgrids_kvec(3, 1, 1, kvecs);
    return pbc;
}

PeriodicBoundaryData make_reduced_pbc()
{
    PeriodicBoundaryData pbc;
    pbc.set_latvec({1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0});
    const std::vector<double> kvecs_ibz{0.0, 0.0, 0.0, TWO_PI / 3.0, 0.0, 0.0};
    const std::vector<std::vector<Vector3_Order<double>>> full_kstars{
        {{0.0, 0.0, 0.0}}, {{1.0 / 3.0, 0.0, 0.0}, {-1.0 / 3.0, 0.0, 0.0}}};
    pbc.set_irreducible_kgrids_kvec(3, 1, 1, kvecs_ibz, full_kstars);
    return pbc;
}

SymmetryContext make_two_atom_inversion_context(const PeriodicBoundaryData& pbc)
{
    SymmetryContext ctx;
    const Matrix3 lattice(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
    ctx.set_crystal_structure(lattice, lattice, {{0, 0}, {1, 0}},
                              {{0, {0.25, 0.0, 0.0}}, {1, {0.75, 0.0, 0.0}}});

    SymmetryOperation identity;
    identity.rotation.Identity();
    identity.translation = {0.0, 0.0, 0.0};

    SymmetryOperation inversion;
    inversion.rotation = Matrix3(-1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
    inversion.translation = {0.0, 0.0, 0.0};

    ctx.set_rspace_operations({identity, inversion});
    ctx.set_available();
    ctx.build_periodic_mappings(pbc, pbc.Rlist);
    ctx.build_rsh_rotations(
        {-1, 0, LIBRPA_ANGULAR_ORDER_NATURAL, LIBRPA_RSH_COEFF_1_M, LIBRPA_RSH_COEFF_1_M}, 0);
    ctx.build_kstar_member_rotations(0);
    return ctx;
}

SymmetryContext make_single_atom_antiunitary_context(const PeriodicBoundaryData& pbc)
{
    SymmetryContext ctx;
    const Matrix3 lattice(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
    ctx.set_crystal_structure(lattice, lattice, {{0, 0}}, {{0, {0.0, 0.0, 0.0}}});

    SymmetryOperation identity;
    identity.rotation.Identity();
    identity.translation = {0.0, 0.0, 0.0};
    ctx.set_rspace_operations({identity});

    SymmetrySpinOperation spin_unitary;
    SymmetrySpinOperation spin_antiunitary;
    spin_antiunitary.antiunitary = true;
    ctx.set_symmetry_spin_operations({spin_unitary, spin_antiunitary}, false);

    ctx.set_available();
    ctx.build_periodic_mappings(pbc, pbc.Rlist);
    ctx.build_rsh_rotations(
        {-1, 0, LIBRPA_ANGULAR_ORDER_NATURAL, LIBRPA_RSH_COEFF_1_M, LIBRPA_RSH_COEFF_1_M}, 0);
    ctx.build_kstar_member_rotations(0);
    return ctx;
}

void assert_full_vr_matches_full(const atpair_R_mat_t& actual, const atpair_R_mat_t& expected,
                                 const double tolerance = 1e-12)
{
    assert(!actual.empty());
    assert(actual.size() == expected.size());
    for (const auto& [atom_i, row] : actual)
    {
        assert(expected.count(atom_i) != 0);
        assert(row.size() == expected.at(atom_i).size());
        for (const auto& [atom_j, R_blocks] : row)
        {
            assert(expected.at(atom_i).count(atom_j) != 0);
            assert(R_blocks.size() == expected.at(atom_i).at(atom_j).size());
            for (const auto& [R, block] : R_blocks)
            {
                assert(expected.at(atom_i).at(atom_j).count(R) != 0);
                const auto& reference = expected.at(atom_i).at(atom_j).at(R);
                assert(block != nullptr && reference != nullptr);
                assert(block->nr == reference->nr && block->nc == reference->nc);
                for (std::size_t index = 0; index != block->size; ++index)
                {
                    const double difference = std::abs(block->c[index] - reference->c[index]);
                    if (difference >= tolerance)
                    {
                        std::cerr << "V(R) mismatch for atom_pair=(" << atom_i << "," << atom_j
                                  << ") R=(" << R.x << "," << R.y << "," << R.z
                                  << ") index=" << index << " actual=" << block->c[index]
                                  << " expected=" << reference->c[index]
                                  << " difference=" << difference << std::endl;
                    }
                    assert(difference < tolerance);
                }
            }
        }
    }
}

void test_vq_to_full_vr_inversion_matches_full_bz()
{
    const auto pbc_full = make_full_pbc();
    const auto pbc_symmetry = make_reduced_pbc();
    auto symmetry_context = make_two_atom_inversion_context(pbc_symmetry);

    AtomicBasis basis_abf(std::vector<std::size_t>{1, 1});
    basis_abf.set_l_shells({{0}, {0}});

    const std::map<atom_t, std::map<atom_t, std::complex<double>>> gamma_blocks{
        {0, {{0, {1.5, 0.0}}, {1, {0.4, 0.0}}}}, {1, {{1, {1.5, 0.0}}}}};
    const std::map<atom_t, std::map<atom_t, std::complex<double>>> representative_blocks{
        {0, {{0, {2.1, 0.0}}, {1, {-0.7, 0.5}}}}, {1, {{1, {1.4, 0.0}}}}};
    const std::map<atom_t, std::map<atom_t, std::complex<double>>> inversion_blocks{
        {0, {{0, {1.4, 0.0}}, {1, {-0.7, -0.5}}}}, {1, {{1, {2.1, 0.0}}}}};

    atpair_k_cplx_mat_t vq_symmetry;
    add_scalar_vq_blocks(vq_symmetry, pbc_symmetry.klist.at(0), gamma_blocks);
    add_scalar_vq_blocks(vq_symmetry, pbc_symmetry.klist.at(1), representative_blocks);
    scale_replicated_collective_input(vq_symmetry);

    atpair_k_cplx_mat_t vq_full;
    add_scalar_vq_blocks(vq_full, pbc_full.klist.at(0), gamma_blocks);
    add_scalar_vq_blocks(vq_full, pbc_full.klist.at(1), representative_blocks);
    add_scalar_vq_blocks(vq_full, pbc_full.klist.at(2), inversion_blocks);

    SymmetryContext no_symmetry;
    bool full_result_replicated = true;
    const auto expected = FT_Vq(global::mpi_comm_global_h, basis_abf, no_symmetry, vq_full,
                                pbc_full, false, false, &full_result_replicated);
    bool symmetry_result_replicated = false;
    const auto actual = FT_Vq(global::mpi_comm_global_h, basis_abf, symmetry_context, vq_symmetry,
                              pbc_symmetry, false, true, &symmetry_result_replicated);

    assert(!full_result_replicated);
    assert(symmetry_result_replicated);
    assert_full_vr_matches_full(actual, expected);
}

void test_vq_to_full_vr_antiunitary_matches_full_bz()
{
    const auto pbc_full = make_full_pbc();
    const auto pbc_symmetry = make_reduced_pbc();
    auto symmetry_context = make_single_atom_antiunitary_context(pbc_symmetry);

    AtomicBasis basis_abf(std::vector<std::size_t>{1});
    basis_abf.set_l_shells({{0}});

    atpair_k_cplx_mat_t vq_symmetry;
    add_scalar_vq_block(vq_symmetry, 0, 0, pbc_symmetry.klist.at(0), {1.5, 0.0});
    add_scalar_vq_block(vq_symmetry, 0, 0, pbc_symmetry.klist.at(1), {2.1, -0.7});
    scale_replicated_collective_input(vq_symmetry);

    atpair_k_cplx_mat_t vq_full;
    add_scalar_vq_block(vq_full, 0, 0, pbc_full.klist.at(0), {1.5, 0.0});
    add_scalar_vq_block(vq_full, 0, 0, pbc_full.klist.at(1), {2.1, -0.7});
    add_scalar_vq_block(vq_full, 0, 0, pbc_full.klist.at(2), {2.1, 0.7});

    SymmetryContext no_symmetry;
    const auto expected =
        FT_Vq(global::mpi_comm_global_h, basis_abf, no_symmetry, vq_full, pbc_full, false, false);
    const auto actual = FT_Vq(global::mpi_comm_global_h, basis_abf, symmetry_context, vq_symmetry,
                              pbc_symmetry, false, true);

    assert_full_vr_matches_full(actual, expected);
}

}  // namespace

int main(int argc, char* argv[])
{
    int provided = 0;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    global::init_global_mpi(MPI_COMM_WORLD);
    global::init_global_io(false, "stdout", false);

    test_vq_to_full_vr_inversion_matches_full_bz();
    test_vq_to_full_vr_antiunitary_matches_full_bz();

    global::finalize_global_io();
    global::finalize_global_mpi();
    MPI_Finalize();
    return 0;
}
