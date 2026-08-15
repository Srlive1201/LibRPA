#include <array>
#include <cassert>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <valarray>

#include "../core/chi0.h"
#include "../core/dielecmodel.h"
#include "../core/epsilon.h"
#include "../core/qpoint_view.h"
#include "../gpu/la_connector.h"
#include "../core/meanfield_mpi.h"
#include "../io/global_io.h"
#include "../math/utils_matrix_m_mpi.h"
#include "../mpi/base_blacs.h"
#include "../mpi/base_mpi.h"
#include "../mpi/kpoint_blacs_parallel_context.h"
#include "../utils/constants.h"
#include "mpi_test_config.h"

using librpa_int::ArrayDesc;
using librpa_int::AtomicBasis;
using librpa_int::atpair_k_cplx_mat_t;
using librpa_int::BlacsCtxtHandler;
using librpa_int::ComplexMatrix;
using librpa_int::diele_func;
using librpa_int::init_local_mat;
using librpa_int::KPointBlacsParallelContext;
using librpa_int::KPointBlacsProcessShape;
using librpa_int::MAJOR;
using librpa_int::Matrix3;
using librpa_int::Matz;
using librpa_int::matrix_m;
using librpa_int::MeanField;
using librpa_int::PeriodicBoundaryData;
using librpa_int::SpeciesBasisLayout;
using librpa_int::SymmetryContext;
using librpa_int::SymmetryKAtomRotation;
using librpa_int::SymmetryKStarMember;
using librpa_int::SymmetryOperation;
using librpa_int::Vector3_Order;
using librpa_int::atom_t;

namespace
{

void test_rspace_symmetry_requires_complete_band_space()
{
    MeanField truncated(1, 1, 2, 3);
    assert(!librpa_int::rspace_symmetry_has_complete_band_space(truncated, -1));
    assert(!librpa_int::rspace_symmetry_has_complete_band_space(truncated, 3));

    MeanField complete(1, 1, 3, 3);
    assert(librpa_int::rspace_symmetry_has_complete_band_space(complete, -1));
    assert(!librpa_int::rspace_symmetry_has_complete_band_space(complete, 2));
}

void assert_complex_close(const std::complex<double> &actual, const std::complex<double> &expected,
                          const double tolerance)
{
    if (std::abs(actual - expected) >= tolerance)
    {
        std::cerr << "actual=" << actual << " expected=" << expected
                  << " diff=" << std::abs(actual - expected) << std::endl;
        assert(false);
    }
}

void require_double_close(const double actual, const double expected, const double tolerance)
{
    if (std::abs(actual - expected) >= tolerance)
    {
        std::cerr << "actual=" << actual << " expected=" << expected
                  << " diff=" << std::abs(actual - expected) << std::endl;
        std::abort();
    }
}

void fill_distributed_matrix(
    matrix_m<std::complex<double>> &matrix, const ArrayDesc &desc,
    const std::vector<std::vector<std::complex<double>>> &values)
{
    assert(static_cast<int>(values.size()) == desc.m());
    for (int i = 0; i != desc.m(); ++i)
    {
        assert(static_cast<int>(values[i].size()) == desc.n());
        const int ilo = desc.indx_g2l_r(i);
        if (ilo < 0) continue;
        for (int j = 0; j != desc.n(); ++j)
        {
            const int jlo = desc.indx_g2l_c(j);
            if (jlo >= 0) matrix(ilo, jlo) = values[i][j];
        }
    }
}

void fill_distributed_matrix(
    matrix_m<std::complex<double>> &matrix, const ArrayDesc &desc,
    const matrix_m<std::complex<double>> &values)
{
    assert(values.nr() == desc.m() && values.nc() == desc.n());
    for (int i_local = 0; i_local != desc.m_loc(); ++i_local)
    {
        const int i_global = desc.indx_l2g_r(i_local);
        for (int j_local = 0; j_local != desc.n_loc(); ++j_local)
        {
            const int j_global = desc.indx_l2g_c(j_local);
            matrix(i_local, j_local) = values(i_global, j_global);
        }
    }
}

void verify_distributed_inverse(
    const std::vector<std::vector<std::complex<double>>> &values,
    const BlacsCtxtHandler &blacs_h, const int block_size, const bool use_cholesky,
    const bool expect_empty_local_rank)
{
    const int n = static_cast<int>(values.size());
    ArrayDesc desc(blacs_h);
    desc.init(n, n, block_size, block_size, 0, 0);

    auto original = init_local_mat<std::complex<double>>(desc, MAJOR::COL);
    fill_distributed_matrix(original, desc, values);
    auto inverse = original.copy();

    int local_is_empty = inverse.size() == 0 ? 1 : 0;
    int any_empty = 0;
    MPI_Allreduce(&local_is_empty, &any_empty, 1, MPI_INT, MPI_MAX, desc.comm());
    if (expect_empty_local_rank) assert(any_empty == 1);

    librpa_int::invert_headwing_body_with_identity_solve(
        inverse, desc, blacs_h, use_cholesky, false);

    auto product = init_local_mat<std::complex<double>>(desc, MAJOR::COL);
    librpa_int::ScalapackConnector::pgemm_f(
        'N', 'N', n, n, n, std::complex<double>{1.0, 0.0}, original.ptr(), 1, 1,
        desc.desc, inverse.ptr(), 1, 1, desc.desc, std::complex<double>{0.0, 0.0},
        product.ptr(), 1, 1, desc.desc);

    double local_max_error = 0.0;
    for (int ilo = 0; ilo != desc.m_loc(); ++ilo)
    {
        const int i = desc.indx_l2g_r(ilo);
        for (int jlo = 0; jlo != desc.n_loc(); ++jlo)
        {
            const int j = desc.indx_l2g_c(jlo);
            const std::complex<double> expected = i == j ? 1.0 : 0.0;
            local_max_error =
                std::max(local_max_error, std::abs(product(ilo, jlo) - expected));
        }
    }
    double max_error = 0.0;
    MPI_Allreduce(&local_max_error, &max_error, 1, MPI_DOUBLE, MPI_MAX, desc.comm());
    require_double_close(max_error, 0.0, 1.0e-10);
}

void test_headwing_body_inverse_uses_identity_solve(const BlacsCtxtHandler &square_blacs_h)
{
    const std::vector<std::vector<std::complex<double>>> pivoting_matrix{
        {{0.0, 0.0}, {1.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}},
        {{1.0, 0.0}, {4.0, 0.0}, {1.0, 0.2}, {0.0, 0.0}},
        {{0.0, 0.0}, {1.0, -0.2}, {5.0, 0.0}, {1.0, 0.0}},
        {{0.0, 0.0}, {0.0, 0.0}, {1.0, 0.0}, {6.0, 0.0}}};
    verify_distributed_inverse(pivoting_matrix, square_blacs_h, 1, false, false);

    BlacsCtxtHandler horizontal_blacs_h(MPI_COMM_WORLD);
    horizontal_blacs_h.init();
    horizontal_blacs_h.set_horizontal_grid();
    const std::vector<std::vector<std::complex<double>>> positive_definite_matrix{
        {{4.0, 0.0}, {1.0, 1.0}},
        {{1.0, -1.0}, {3.0, 0.0}}};
    verify_distributed_inverse(positive_definite_matrix, horizontal_blacs_h, 2, true,
                               horizontal_blacs_h.nprocs > 1);
}

void test_gamma_head_rank_one_matches_coulomb_basis_overwrite(
    const BlacsCtxtHandler &blacs_h)
{
    constexpr int n = 4;
    const double half = 0.5;
    const matrix_m<std::complex<double>> eigenvectors(
        {{{half, 0.0}, {half, 0.0}, {half, 0.0}, {half, 0.0}},
         {{half, 0.0}, {-half, 0.0}, {half, 0.0}, {-half, 0.0}},
         {{half, 0.0}, {half, 0.0}, {-half, 0.0}, {-half, 0.0}},
         {{half, 0.0}, {-half, 0.0}, {-half, 0.0}, {half, 0.0}}},
        MAJOR::COL);
    const std::array<double, n> eigenvalues{{9.0, 4.0, 1.0, 0.0}};
    const matrix_m<std::complex<double>> response(
        {{{0.13, 0.0}, {0.02, 0.03}, {-0.01, 0.02}, {0.04, -0.01}},
         {{0.02, -0.03}, {0.18, 0.0}, {0.03, -0.02}, {-0.02, 0.01}},
         {{-0.01, -0.02}, {0.03, 0.02}, {0.16, 0.0}, {0.01, 0.04}},
         {{0.04, 0.01}, {-0.02, -0.01}, {0.01, -0.04}, {0.11, 0.0}}},
        MAJOR::COL);
    const std::complex<double> corrected_head{2.35, -0.015};

    const auto eigenvectors_h = eigenvectors.get_transpose(true);
    auto scaled_eigenvectors = eigenvectors.copy();
    for (int i = 0; i != n; ++i)
        scaled_eigenvectors.scale_col(i, std::sqrt(eigenvalues.at(i)));
    const auto sqrt_coulomb = scaled_eigenvectors * eigenvectors_h;

    auto epsilon_direct = sqrt_coulomb * response * sqrt_coulomb;
    epsilon_direct *= -1.0;
    for (int i = 0; i != n; ++i) epsilon_direct(i, i) += 1.0;

    const auto epsilon_eigenbasis = eigenvectors_h * epsilon_direct * eigenvectors;
    const auto uncorrected_head_reference = epsilon_eigenbasis(0, 0);
    auto response_eigenbasis = eigenvectors_h * response * eigenvectors;
    for (int i = 0; i != n; ++i)
    {
        for (int j = 0; j != n; ++j)
        {
            response_eigenbasis(i, j) *=
                -std::sqrt(eigenvalues.at(i) * eigenvalues.at(j));
        }
    }
    response_eigenbasis(0, 0) = corrected_head - 1.0;
    auto epsilon_reference = eigenvectors * response_eigenbasis * eigenvectors_h;
    for (int i = 0; i != n; ++i) epsilon_reference(i, i) += 1.0;

    ArrayDesc desc(blacs_h);
    desc.init(n, n, 1, 2, 0, 0);
    auto eigen_local = init_local_mat<std::complex<double>>(desc, MAJOR::COL);
    auto epsilon_local = init_local_mat<std::complex<double>>(desc, MAJOR::COL);
    auto scratch_local = init_local_mat<std::complex<double>>(desc, MAJOR::COL);
    fill_distributed_matrix(eigen_local, desc, eigenvectors);
    fill_distributed_matrix(epsilon_local, desc, epsilon_direct);

    ArrayDesc desc_1x1(blacs_h);
    desc_1x1.init(1, 1, 2, 2, 0, 0);
    auto h_local = init_local_mat<std::complex<double>>(desc_1x1, MAJOR::COL);
    if (desc_1x1.m_loc() == 0 || desc_1x1.n_loc() == 0) h_local.resize(1, 1);
    if (desc_1x1.m_loc() != 0 && desc_1x1.n_loc() != 0) h_local(0, 0) = {0.0, 0.0};

    // y = epsilon0 * x1
    librpa_int::LaConnector::pgemm(
        'N', 'N', n, 1, n, std::complex<double>{1.0, 0.0},
        epsilon_local.ptr(), 1, 1, desc,
        eigen_local.ptr(), 1, 1, desc,
        std::complex<double>{0.0, 0.0},
        scratch_local.ptr(), 1, 1, desc);

    // h = x1^H * y
    librpa_int::LaConnector::pgemm(
        'C', 'N', 1, 1, n, std::complex<double>{1.0, 0.0},
        eigen_local.ptr(), 1, 1, desc,
        scratch_local.ptr(), 1, 1, desc,
        std::complex<double>{0.0, 0.0},
        h_local.ptr(), 1, 1, desc_1x1);

    std::complex<double> h_scalar{0.0, 0.0};
    if (desc_1x1.is_src()) h_scalar = h_local(0, 0);
    const int h_root = blacs_h.get_pnum(0, 0);
    MPI_Bcast(&h_scalar, 1, MPI_CXX_DOUBLE_COMPLEX, h_root, desc_1x1.comm());

    assert_complex_close(h_scalar, uncorrected_head_reference, 1.0e-12);

    // epsilon += (H - h) * x1 * x1^H
    const std::complex<double> coeff = corrected_head - h_scalar;
    librpa_int::LaConnector::pgemm(
        'N', 'C', n, n, 1, coeff,
        eigen_local.ptr(), 1, 1, desc,
        eigen_local.ptr(), 1, 1, desc,
        std::complex<double>{1.0, 0.0},
        epsilon_local.ptr(), 1, 1, desc);

    double local_max_error = 0.0;
    for (int i_local = 0; i_local != desc.m_loc(); ++i_local)
    {
        const int i_global = desc.indx_l2g_r(i_local);
        for (int j_local = 0; j_local != desc.n_loc(); ++j_local)
        {
            const int j_global = desc.indx_l2g_c(j_local);
            local_max_error = std::max(
                local_max_error,
                std::abs(epsilon_local(i_local, j_local)
                         - epsilon_reference(i_global, j_global)));
        }
    }
    double max_error = 0.0;
    MPI_Allreduce(&local_max_error, &max_error, 1, MPI_DOUBLE, MPI_MAX, desc.comm());
    require_double_close(max_error, 0.0, 1.0e-11);

    // Verify that modifying eigenvector columns other than column zero
    // cannot affect h or the rank-one update.
    auto eigen_modified = eigen_local.copy();
    for (int col = 1; col != n; ++col)
    {
        const int col_local = desc.indx_g2l_c(col);
        if (col_local < 0) continue;
        for (int i_local = 0; i_local != desc.m_loc(); ++i_local)
            eigen_modified(i_local, col_local) *= std::complex<double>{0.0, 2.0};
    }
    auto epsilon_modified = init_local_mat<std::complex<double>>(desc, MAJOR::COL);
    fill_distributed_matrix(epsilon_modified, desc, epsilon_direct);
    auto scratch2 = init_local_mat<std::complex<double>>(desc, MAJOR::COL);
    auto h2_local = init_local_mat<std::complex<double>>(desc_1x1, MAJOR::COL);
    if (desc_1x1.m_loc() == 0 || desc_1x1.n_loc() == 0) h2_local.resize(1, 1);
    if (desc_1x1.m_loc() != 0 && desc_1x1.n_loc() != 0) h2_local(0, 0) = {0.0, 0.0};

    librpa_int::LaConnector::pgemm(
        'N', 'N', n, 1, n, std::complex<double>{1.0, 0.0},
        epsilon_modified.ptr(), 1, 1, desc,
        eigen_modified.ptr(), 1, 1, desc,
        std::complex<double>{0.0, 0.0},
        scratch2.ptr(), 1, 1, desc);
    librpa_int::LaConnector::pgemm(
        'C', 'N', 1, 1, n, std::complex<double>{1.0, 0.0},
        eigen_modified.ptr(), 1, 1, desc,
        scratch2.ptr(), 1, 1, desc,
        std::complex<double>{0.0, 0.0},
        h2_local.ptr(), 1, 1, desc_1x1);

    std::complex<double> h2_scalar{0.0, 0.0};
    if (desc_1x1.is_src()) h2_scalar = h2_local(0, 0);
    MPI_Bcast(&h2_scalar, 1, MPI_CXX_DOUBLE_COMPLEX, h_root, desc_1x1.comm());
    assert_complex_close(h2_scalar, h_scalar, 1.0e-12);

    librpa_int::LaConnector::pgemm(
        'N', 'C', n, n, 1, coeff,
        eigen_modified.ptr(), 1, 1, desc,
        eigen_modified.ptr(), 1, 1, desc,
        std::complex<double>{1.0, 0.0},
        epsilon_modified.ptr(), 1, 1, desc);

    double local_mod_error = 0.0;
    for (int i_local = 0; i_local != desc.m_loc(); ++i_local)
    {
        for (int j_local = 0; j_local != desc.n_loc(); ++j_local)
        {
            local_mod_error = std::max(
                local_mod_error,
                std::abs(epsilon_modified(i_local, j_local)
                         - epsilon_local(i_local, j_local)));
        }
    }
    double mod_error = 0.0;
    MPI_Allreduce(&local_mod_error, &mod_error, 1, MPI_DOUBLE, MPI_MAX, desc.comm());
    require_double_close(mod_error, 0.0, 1.0e-11);
}

void test_gamma_head_rank_one_handles_empty_local_blocks(const BlacsCtxtHandler &blacs_h)
{
    ArrayDesc desc(blacs_h);
    desc.init(1, 1, 1, 1, 0, 0);
    auto eigen_local = init_local_mat<std::complex<double>>(desc, MAJOR::COL);
    auto epsilon_local = init_local_mat<std::complex<double>>(desc, MAJOR::COL);
    auto scratch_local = init_local_mat<std::complex<double>>(desc, MAJOR::COL);
    if (desc.m_loc() == 0 || desc.n_loc() == 0)
    {
        eigen_local.resize(1, 1);
        epsilon_local.resize(1, 1);
        scratch_local.resize(1, 1);
    }
    if (desc.m_loc() != 0 && desc.n_loc() != 0)
    {
        eigen_local(0, 0) = 1.0;
        epsilon_local(0, 0) = 0.8;
    }

    ArrayDesc desc_1x1(blacs_h);
    desc_1x1.init(1, 1, 1, 1, 0, 0);
    auto h_local = init_local_mat<std::complex<double>>(desc_1x1, MAJOR::COL);
    if (desc_1x1.m_loc() == 0 || desc_1x1.n_loc() == 0) h_local.resize(1, 1);
    if (desc_1x1.m_loc() != 0 && desc_1x1.n_loc() != 0) h_local(0, 0) = {0.0, 0.0};

    // y = epsilon0 * x1
    librpa_int::LaConnector::pgemm(
        'N', 'N', 1, 1, 1, std::complex<double>{1.0, 0.0},
        epsilon_local.ptr(), 1, 1, desc,
        eigen_local.ptr(), 1, 1, desc,
        std::complex<double>{0.0, 0.0},
        scratch_local.ptr(), 1, 1, desc);

    // h = x1^H * y
    librpa_int::LaConnector::pgemm(
        'C', 'N', 1, 1, 1, std::complex<double>{1.0, 0.0},
        eigen_local.ptr(), 1, 1, desc,
        scratch_local.ptr(), 1, 1, desc,
        std::complex<double>{0.0, 0.0},
        h_local.ptr(), 1, 1, desc_1x1);

    std::complex<double> h_scalar{0.0, 0.0};
    if (desc_1x1.is_src()) h_scalar = h_local(0, 0);
    const int h_root = blacs_h.get_pnum(0, 0);
    MPI_Bcast(&h_scalar, 1, MPI_CXX_DOUBLE_COMPLEX, h_root, desc_1x1.comm());

    assert_complex_close(h_scalar, std::complex<double>{0.8, 0.0}, 1.0e-12);

    // epsilon += (2.1 - h) * x1 * x1^H
    const std::complex<double> coeff = std::complex<double>{2.1, 0.0} - h_scalar;
    librpa_int::LaConnector::pgemm(
        'N', 'C', 1, 1, 1, coeff,
        eigen_local.ptr(), 1, 1, desc,
        eigen_local.ptr(), 1, 1, desc,
        std::complex<double>{1.0, 0.0},
        epsilon_local.ptr(), 1, 1, desc);

    int local_is_empty = (desc.m_loc() == 0 || desc.n_loc() == 0) ? 1 : 0;
    int any_empty = 0;
    MPI_Allreduce(&local_is_empty, &any_empty, 1, MPI_INT, MPI_MAX, desc.comm());
    if (desc.nprocs() > 1) assert(any_empty == 1);
    if (desc.m_loc() != 0 && desc.n_loc() != 0)
        assert_complex_close(epsilon_local(0, 0), std::complex<double>{2.1, 0.0}, 1.0e-12);
}

RI::Tensor<double> make_scalar_cs_tensor(const double value)
{
    auto data = std::make_shared<std::valarray<double>>(1);
    (*data)[0] = value;
    return RI::Tensor<double>({1UL, 1UL, 1UL}, data);
}

void test_kpoint_coordinate_mapping_selects_active_klist_from_full_source()
{
    const std::vector<Vector3_Order<double>> pyatb_full_kpoints{
        {0.0, 0.0, 0.0}, {0.125, 0.0, 0.0}, {0.25, 0.0, 0.0}, {0.375, 0.0, 0.0},
        {0.0, 0.125, 0.0}, {0.125, 0.125, 0.0}, {0.875, 0.0, 0.0}};
    const std::vector<Vector3_Order<double>> active_kpoints{
        {0.0, 0.0, 0.0}, {0.25, 0.0, 0.0}, {0.125, 0.125, 0.0}};

    const auto mapping = librpa_int::map_kpoints_by_coordinates(active_kpoints, pyatb_full_kpoints);

    assert((mapping == std::vector<int>{0, 2, 5}));

    const std::vector<Vector3_Order<double>> wrapped_active_kpoints{{-0.125, 0.0, 0.0}};
    const auto wrapped_mapping =
        librpa_int::map_kpoints_by_coordinates(wrapped_active_kpoints, pyatb_full_kpoints);
    assert((wrapped_mapping == std::vector<int>{6}));
}

void test_kstar_velocity_mapping_preserves_member_order_and_periodic_gauge()
{
    SymmetryContext ctx;
    librpa_int::SymmetryKStar gamma;
    gamma.k_ibz = {0.0, 0.0, 0.0};
    librpa_int::SymmetryKStarMember gamma_first;
    gamma_first.k_bz = {0.0, 0.0, 0.0};
    librpa_int::SymmetryKStarMember gamma_second;
    gamma_second.k_bz = {0.5, 0.0, 0.0};
    gamma.members = {gamma_first, gamma_second};
    librpa_int::SymmetryKStar quarter;
    quarter.k_ibz = {0.25, 0.0, 0.0};
    librpa_int::SymmetryKStarMember quarter_first;
    quarter_first.k_bz = {0.25, 0.0, 0.0};
    librpa_int::SymmetryKStarMember quarter_second;
    quarter_second.k_bz = {-0.25, 0.0, 0.0};
    quarter.members = {quarter_first, quarter_second};
    ctx.kstars = {gamma, quarter};

    const std::vector<Vector3_Order<double>> ibz_kpoints{{0.25, 0.0, 0.0}, {0.0, 0.0, 0.0}};
    const std::vector<Vector3_Order<double>> full_bz_kpoints{
        {0.5, 0.0, 0.0}, {0.75, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.25, 0.0, 0.0}};

    const auto mapping = librpa_int::map_symmetry_kstar_members_to_source_kpoints(
        ctx, ibz_kpoints, full_bz_kpoints);

    assert((mapping == std::vector<std::vector<int>>{{3, 1}, {2, 0}}));
}

void test_replace_rpa_response_headwing_replaces_only_singular_channels(
    const BlacsCtxtHandler &blacs_h)
{
    ArrayDesc desc(blacs_h);
    desc.init_square_blk(4, 4, 0, 0);

    auto response = init_local_mat<std::complex<double>>(desc, MAJOR::COL);
    for (int i = 0; i != 4; ++i)
    {
        const int ilo = desc.indx_g2l_r(i);
        if (ilo < 0) continue;
        for (int j = 0; j != 4; ++j)
        {
            const int jlo = desc.indx_g2l_c(j);
            if (jlo < 0) continue;
            response(ilo, jlo) =
                std::complex<double>(0.01 * (i + 1) + 0.02 * (j + 1), 0.001 * (i - j));
        }
    }

    const matrix_m<std::complex<double>> head(
        std::vector<std::vector<std::complex<double>>>{
            {std::complex<double>{2.0, 0.0}, std::complex<double>{0.2, 0.1},
             std::complex<double>{0.3, -0.1}},
            {std::complex<double>{0.2, -0.1}, std::complex<double>{2.2, 0.0},
             std::complex<double>{0.4, 0.2}},
            {std::complex<double>{0.3, 0.1}, std::complex<double>{0.4, -0.2},
             std::complex<double>{2.4, 0.0}}},
        MAJOR::COL);
    const matrix_m<std::complex<double>> wing(
        std::vector<std::vector<std::complex<double>>>{
            {std::complex<double>{0.11, 0.01}, std::complex<double>{0.12, 0.02},
             std::complex<double>{0.13, 0.03}},
            {std::complex<double>{0.21, 0.04}, std::complex<double>{0.22, 0.05},
             std::complex<double>{0.23, 0.06}},
            {std::complex<double>{0.31, 0.07}, std::complex<double>{0.32, 0.08},
             std::complex<double>{0.33, 0.09}}},
        MAJOR::COL);

    librpa_int::replace_rpa_response_headwing(response, head, wing, desc);

    const int head_row = desc.indx_g2l_r(0);
    const int head_col = desc.indx_g2l_c(0);
    if (head_row >= 0 && head_col >= 0)
    {
        assert_complex_close(response(head_row, head_col), std::complex<double>{2.2, 0.0}, 1e-12);
    }

    for (int lambda = 1; lambda != 4; ++lambda)
    {
        std::complex<double> expected_wing = 0.0;
        for (int alpha = 0; alpha != 3; ++alpha)
        {
            expected_wing += wing(lambda - 1, alpha);
        }
        expected_wing /= 3.0;

        const int row_body = desc.indx_g2l_r(lambda);
        const int col_head = desc.indx_g2l_c(0);
        if (row_body >= 0 && col_head >= 0)
        {
            assert_complex_close(response(row_body, col_head), expected_wing, 1e-12);
        }

        const int row_head = desc.indx_g2l_r(0);
        const int col_body = desc.indx_g2l_c(lambda);
        if (row_head >= 0 && col_body >= 0)
        {
            assert_complex_close(response(row_head, col_body), std::conj(expected_wing), 1e-12);
        }
    }

    const int body_row = desc.indx_g2l_r(2);
    const int body_col = desc.indx_g2l_c(3);
    if (body_row >= 0 && body_col >= 0)
    {
        assert_complex_close(response(body_row, body_col), std::complex<double>{0.11, -0.001},
                             1e-12);
    }
}

void test_replace_rpa_response_head_only_keeps_numeric_wings(const BlacsCtxtHandler &blacs_h)
{
    ArrayDesc desc(blacs_h);
    desc.init_square_blk(4, 4, 0, 0);

    auto response = init_local_mat<std::complex<double>>(desc, MAJOR::COL);
    matrix_m<std::complex<double>> original(4, 4, MAJOR::COL);
    for (int i = 0; i != 4; ++i)
    {
        const int ilo = desc.indx_g2l_r(i);
        for (int j = 0; j != 4; ++j)
        {
            const auto value =
                std::complex<double>(0.1 * (i + 1) + 0.01 * (j + 1), 0.001 * (i - j));
            original(i, j) = value;
            if (ilo < 0) continue;
            const int jlo = desc.indx_g2l_c(j);
            if (jlo < 0) continue;
            response(ilo, jlo) = value;
        }
    }

    const matrix_m<std::complex<double>> chi0v_head(
        std::vector<std::vector<std::complex<double>>>{
            {std::complex<double>{0.21, 0.0}, std::complex<double>{0.01, 0.02},
             std::complex<double>{-0.03, 0.04}},
            {std::complex<double>{0.05, -0.01}, std::complex<double>{0.24, 0.0},
             std::complex<double>{0.07, 0.03}},
            {std::complex<double>{-0.02, -0.04}, std::complex<double>{0.08, -0.03},
             std::complex<double>{0.27, 0.0}}},
        MAJOR::COL);

    librpa_int::replace_rpa_response_head_only(response, chi0v_head, desc);

    const auto expected_head = (chi0v_head(0, 0) + chi0v_head(1, 1) + chi0v_head(2, 2)) / 3.0;
    for (int i = 0; i != 4; ++i)
    {
        const int ilo = desc.indx_g2l_r(i);
        if (ilo < 0) continue;
        for (int j = 0; j != 4; ++j)
        {
            const int jlo = desc.indx_g2l_c(j);
            if (jlo < 0) continue;
            const auto expected = (i == 0 && j == 0) ? expected_head : original(i, j);
            assert_complex_close(response(ilo, jlo), expected, 1e-12);
        }
    }
}

void test_rpa_trace_log_average_uses_directional_head_and_wing()
{
    const matrix_m<std::complex<double>> head(
        std::vector<std::vector<std::complex<double>>>{
            {std::complex<double>{0.2, 0.0}, std::complex<double>{0.0, 0.0},
             std::complex<double>{0.0, 0.0}},
            {std::complex<double>{0.0, 0.0}, std::complex<double>{0.3, 0.0},
             std::complex<double>{0.0, 0.0}},
            {std::complex<double>{0.0, 0.0}, std::complex<double>{0.0, 0.0},
             std::complex<double>{0.4, 0.0}}},
        MAJOR::COL);
    const std::complex<double> body{0.1, 0.0};
    const std::array<std::complex<double>, 3> wing{std::complex<double>{0.05, 0.0},
                                                   std::complex<double>{0.02, 0.0},
                                                   std::complex<double>{0.01, 0.0}};
    const std::complex<double> body_inv = 1.0 / (1.0 - body);
    const matrix_m<std::complex<double>> schur_l(
        std::vector<std::vector<std::complex<double>>>{
            {1.0 - head(0, 0) - std::conj(wing[0]) * body_inv * wing[0],
             -std::conj(wing[0]) * body_inv * wing[1], -std::conj(wing[0]) * body_inv * wing[2]},
            {-std::conj(wing[1]) * body_inv * wing[0],
             1.0 - head(1, 1) - std::conj(wing[1]) * body_inv * wing[1],
             -std::conj(wing[1]) * body_inv * wing[2]},
            {-std::conj(wing[2]) * body_inv * wing[0], -std::conj(wing[2]) * body_inv * wing[1],
             1.0 - head(2, 2) - std::conj(wing[2]) * body_inv * wing[2]}},
        MAJOR::COL);
    const std::vector<double> qx{1.0, 0.0};
    const std::vector<double> qy{0.0, 1.0};
    const std::vector<double> qz{0.0, 0.0};
    const std::vector<double> weights{0.2, 0.6};
    const std::complex<double> trace_body = body;
    const std::complex<double> logdet_body = std::log(1.0 - body);

    const auto actual = librpa_int::compute_rpa_chi0v_headwing_trace_log_average(
        head, schur_l, trace_body, logdet_body, qx, qy, qz, weights);
    const auto direct_trace_log = [&](const double nx, const double ny, const double nz)
    {
        const auto directional_head = nx * (nx * head(0, 0) + ny * head(0, 1) + nz * head(0, 2)) +
                                      ny * (nx * head(1, 0) + ny * head(1, 1) + nz * head(1, 2)) +
                                      nz * (nx * head(2, 0) + ny * head(2, 1) + nz * head(2, 2));
        const auto directional_wing = nx * wing[0] + ny * wing[1] + nz * wing[2];
        const auto direct_det = (1.0 - directional_head) * (1.0 - body) -
                                std::conj(directional_wing) * directional_wing;
        return directional_head + body + std::log(direct_det);
    };
    const auto expected = weights[0] * direct_trace_log(qx[0], qy[0], qz[0]) +
                          weights[1] * direct_trace_log(qx[1], qy[1], qz[1]);

    assert_complex_close(actual, expected, 1e-12);
}

void test_rpa_headwing_regular_body_start_channel()
{
    librpa_int::RpaHeadwingSettings settings;

    settings.use_2d_dielectric = false;
    settings.rpa_headwing_body_start = 0;
    assert(librpa_int::rpa_headwing_regular_body_start_channel(settings) == 1);

    settings.use_2d_dielectric = true;
    settings.rpa_headwing_body_start = 0;
    assert(librpa_int::rpa_headwing_regular_body_start_channel(settings) == 1);

    settings.use_2d_dielectric = false;
    settings.rpa_headwing_body_start = 1;
    assert(librpa_int::rpa_headwing_regular_body_start_channel(settings) == 1);

    settings.use_2d_dielectric = true;
    settings.rpa_headwing_body_start = 4;
    assert(librpa_int::rpa_headwing_regular_body_start_channel(settings) == 4);
}

void test_rpa_headwing_gamma_cell_volume_uses_reciprocal_lattice()
{
    librpa_int::PeriodicBoundaryData pbc;
    pbc.latvec = librpa_int::Matrix3(2.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 5.0);
    pbc.G = librpa_int::Matrix3(0.5, 0.0, 0.0, 0.0, 1.0 / 3.0, 0.0, 0.0, 0.0, 0.2);

    const double vol_3d = librpa_int::rpa_headwing_reciprocal_cell_volume(pbc, false);
    require_double_close(vol_3d, std::abs(pbc.G.Det()), 1e-14);

    const double vol_2d = librpa_int::rpa_headwing_reciprocal_cell_volume(pbc, true);
    const double expected_2d = std::abs(pbc.G.e11 * pbc.G.e22 - pbc.G.e12 * pbc.G.e21);
    require_double_close(vol_2d, expected_2d, 1e-14);

    pbc.set_period(4, 4, 4);
    require_double_close(librpa_int::rpa_headwing_gamma_cell_volume(pbc, false),
                         vol_3d / 64.0, 1e-14);
    require_double_close(librpa_int::rpa_headwing_gamma_cell_volume(pbc, true),
                         vol_2d / 64.0, 1e-14);
}

void test_rpa_chi0v_wing_desc_uses_global_rows(const BlacsCtxtHandler &blacs_h)
{
    ArrayDesc desc_body(blacs_h);
    desc_body.init_square_blk(10, 10, 0, 0);

    ArrayDesc desc_full_wing(blacs_h);
    desc_full_wing.init(11, 3, desc_body.mb(), 1, 0, 0);

    const auto desc_wing = librpa_int::make_rpa_chi0v_wing_desc(
        desc_body, 1, desc_full_wing.m_loc(), desc_full_wing.n_loc());

    if (desc_body.nprows() > 1)
    {
        assert(desc_full_wing.m_loc() < desc_full_wing.m());
    }
    assert(desc_wing.m() == 11);
    assert(desc_wing.n() == 3);
    assert(desc_wing.mb() == desc_body.mb());
    assert(desc_wing.nb() == 1);
    assert(desc_wing.m_loc() == desc_full_wing.m_loc());
    assert(desc_wing.n_loc() == desc_full_wing.n_loc());
}

void test_headwing_spin_weights()
{
    assert(std::abs(librpa_int::headwing_transition_weight(1.0, 0.25, 2, false) - 0.75) < 1e-12);
    assert(std::abs(librpa_int::headwing_spin_prefactor(2, false) - 1.0) < 1e-12);

    assert(std::abs(librpa_int::headwing_transition_weight(1.0, 0.25, 1, true) - 0.75) < 1e-12);
    assert(std::abs(librpa_int::headwing_spin_prefactor(1, true) - 1.0) < 1e-12);

    assert(std::abs(librpa_int::headwing_transition_weight(1.0, 0.25, 1, false) - 0.375) < 1e-12);
    assert(std::abs(librpa_int::headwing_spin_prefactor(1, false) - 2.0) < 1e-12);
}

void test_wing_cartesian_gram_is_invariant_under_row_phases()
{
    ComplexMatrix wing(2, 3);
    wing(0, 0) = {1.0, 1.0};
    wing(0, 1) = {2.0, -1.0};
    wing(0, 2) = {0.0, -1.0};
    wing(1, 0) = {0.5, -2.0};
    wing(1, 1) = {-1.0, 0.25};
    wing(1, 2) = {3.0, 0.5};

    ComplexMatrix phased = wing;
    for (int alpha = 0; alpha != 3; ++alpha)
    {
        phased(0, alpha) *= std::complex<double>{0.0, 1.0};
        phased(1, alpha) *= std::complex<double>{-1.0, 0.0};
    }

    const auto gram = librpa_int::compute_wing_cartesian_gram(wing);
    const auto phased_gram = librpa_int::compute_wing_cartesian_gram(phased);
    for (int alpha = 0; alpha != 3; ++alpha)
    {
        for (int beta = 0; beta != 3; ++beta)
        {
            assert_complex_close(gram.at(alpha).at(beta), phased_gram.at(alpha).at(beta), 1e-12);
        }
    }
    assert_complex_close(gram.at(0).at(0), {6.25, 0.0}, 1e-12);
}

void test_velocity_matrix_initialization()
{
    librpa_int::velocity_matrix_t velocity;
    librpa_int::initialize_velocity_matrix(velocity, 2, 3, 4);

    assert(velocity.size() == 2);
    for (int ispin = 0; ispin != 2; ++ispin)
    {
        assert(velocity[ispin].size() == 3);
        for (int ik = 0; ik != 3; ++ik)
        {
            assert(velocity[ispin][ik].size() == 3);
            for (int alpha = 0; alpha != 3; ++alpha)
            {
                assert(velocity[ispin][ik][alpha].nr == 4);
                assert(velocity[ispin][ik][alpha].nc == 4);
            }
        }
    }
}

void test_headwing_local_kpoints_prefers_kpoint_blacs_context()
{
    const auto all_k = librpa_int::headwing_local_kpoints(4, nullptr);
    assert((all_k == std::vector<int>{0, 1, 2, 3}));

    if (librpa_int::get_mpi_size(MPI_COMM_WORLD) != 4) return;

    KPointBlacsProcessShape shape(2, 2, true);
    KPointBlacsParallelContext kctx(shape, MPI_COMM_WORLD, 4);
    const auto local_k = librpa_int::headwing_local_kpoints(4, &kctx);

    if (kctx.kpoint_group_id() == 0)
        assert((local_k == std::vector<int>{0, 2}));
    else
        assert((local_k == std::vector<int>{1, 3}));

    const auto mismatched = librpa_int::headwing_local_kpoints(5, &kctx);
    assert((mismatched == std::vector<int>{0, 1, 2, 3, 4}));
}

void test_headwing_world_fourier_uses_all_R_blocks_at_nonzero_k()
{
    AtomicBasis basis_wfc({1});
    AtomicBasis basis_abf({1});
    librpa_int::Cs_LRI Cs_data;
    Cs_data.use_libri = true;
    Cs_data.data_libri[0][{0, {0, 0, 0}}] = make_scalar_cs_tensor(2.0);
    Cs_data.data_libri[0][{0, {1, 0, 0}}] = make_scalar_cs_tensor(3.0);
    Cs_data.data_libri[0][{0, {2, 0, 0}}] = make_scalar_cs_tensor(-1.0);

    const auto targets = librpa_int::build_headwing_full_bz_fourier_targets(
        {{0.0, 0.0, 0.0}, {0.25, 0.0, 0.0}});
    const auto Cs_IJ_k =
        librpa_int::fourier_headwing_cs_to_ijk(Cs_data, basis_wfc, basis_abf, targets);

    const auto &blocks = Cs_IJ_k.at(0);
    assert_complex_close(blocks.at({0, 0})(0, 0, 0), {4.0, 0.0}, 1e-12);
    assert_complex_close(blocks.at({0, 1})(0, 0, 0), {3.0, 3.0}, 1e-12);
}

void test_headwing_symmetry_fourier_target_ids_are_deterministic()
{
    librpa_int::SymmetryContext ctx;
    ctx.set_available();

    librpa_int::SymmetryKStar first;
    first.star_index = 0;
    first.k_ibz = {0.0, 0.0, 0.0};
    first.members.resize(2);
    first.members[0].k_bz = {0.0, 0.0, 0.0};
    first.members[1].k_bz = {0.5, 0.0, 0.0};
    ctx.kstars.push_back(first);

    librpa_int::SymmetryKStar second;
    second.star_index = 1;
    second.k_ibz = {0.25, 0.0, 0.0};
    second.members.resize(3);
    second.members[0].k_bz = {0.25, 0.0, 0.0};
    second.members[1].k_bz = {0.0, 0.25, 0.0};
    second.members[2].k_bz = {0.0, 0.0, 0.25};
    ctx.kstars.push_back(second);

    PeriodicBoundaryData pbc;
    const auto flattened = librpa_int::build_headwing_symmetry_fourier_targets(
        ctx, pbc, {first.k_ibz, second.k_ibz});

    assert((flattened.target_ids_by_ibz_member[0] == std::vector<int>{0, 1}));
    assert((flattened.target_ids_by_ibz_member[1] == std::vector<int>{2, 3, 4}));
    assert(flattened.targets.size() == 5);
    for (int target_id = 0; target_id != 5; ++target_id)
        assert(flattened.targets[target_id].target_id == target_id);
    assert(flattened.targets[0].owner_ik == 0);
    assert(flattened.targets[1].owner_ik == 0);
    assert(flattened.targets[2].owner_ik == 1);
    assert(flattened.targets[4].kfrac == second.members[2].k_bz);
}

void test_headwing_ijk_redistribution_is_owner_group_local()
{
    if (librpa_int::get_mpi_size(MPI_COMM_WORLD) != 4) return;

    KPointBlacsProcessShape shape(2, 2, true);
    KPointBlacsParallelContext kctx(shape, MPI_COMM_WORLD, 4);
    const auto desc_nao_nao = kctx.create_array_desc(2, 2, 1, 1);
    const AtomicBasis basis_wfc(std::vector<std::size_t>{1, 1});
    const AtomicBasis basis_abf(std::vector<std::size_t>{1, 1});
    const auto targets = librpa_int::build_headwing_full_bz_fourier_targets(
        {{0.0, 0.0, 0.0}, {0.25, 0.0, 0.0}, {0.5, 0.0, 0.0}, {0.75, 0.0, 0.0}});
    const std::set<int> local_iks(kctx.kpoints_local().begin(), kctx.kpoints_local().end());
    const auto requests = librpa_int::build_headwing_cs_ijk_requests(
        basis_wfc, targets, kctx.kpoints_local(), desc_nao_nao);

    for (const auto &[J, target_id] : requests.second)
    {
        (void)J;
        assert(local_iks.count(targets.at(target_id).owner_ik) == 1);
    }

#ifdef LIBRPA_USE_LIBRI
    librpa_int::Cs_LRI Cs_data;
    Cs_data.use_libri = true;
    if (librpa_int::get_mpi_rank(MPI_COMM_WORLD) == 0)
    {
        for (int I = 0; I != 2; ++I)
            for (int J = 0; J != 2; ++J)
                Cs_data.data_libri[I][{J, {0, 0, 0}}] =
                    make_scalar_cs_tensor(1.0 + 2.0 * I + J);
    }
    const auto Cs_IJ_k = librpa_int::redistribute_headwing_cs_ijk(
        Cs_data, basis_wfc, basis_abf, targets, kctx.kpoints_local(), desc_nao_nao,
        librpa_int::global::mpi_comm_global_h);
    for (const auto &[I, Jtargets] : Cs_IJ_k)
    {
        assert(requests.first.count(I) == 1);
        for (const auto &[Jtarget, tensor] : Jtargets)
        {
            assert(requests.second.count(Jtarget) == 1);
            assert(local_iks.count(targets.at(Jtarget.second).owner_ik) == 1);
            assert(std::abs(tensor(0, 0, 0)) > 0.0);
        }
    }
#endif
}

void test_accumulate_wing_mu_for_pair_matches_original_formula()
{
    const std::vector<double> omega{0.5, 1.25};
    const std::array<std::complex<double>, 3> velocity{std::complex<double>{0.2, -0.1},
                                                       std::complex<double>{-0.3, 0.4},
                                                       std::complex<double>{0.15, 0.05}};
    const std::complex<double> c_mn{0.7, -0.2};
    const double egap = 1.8;
    const double factor1 = 0.6;
    const double factor2 = 0.125;

    std::array<std::complex<double>, 6> accumulated{};
    librpa_int::accumulate_wing_mu_for_pair(omega, velocity, c_mn, egap, factor1, factor2,
                                            accumulated.data());

    for (std::size_t iomega = 0; iomega != omega.size(); ++iomega)
    {
        for (int alpha = 0; alpha != 3; ++alpha)
        {
            const auto denom = omega[iomega] * omega[iomega] + egap * egap;
            const auto expected = factor1 * std::conj(c_mn * velocity[alpha]) / denom +
                                  factor2 * c_mn * velocity[alpha] / denom;
            assert_complex_close(accumulated[iomega * 3 + alpha], expected, 1e-12);
        }
    }
}

SymmetryKStarMember make_headwing_wfc_atom_swap_member(
    const std::complex<double> &rot_0, const std::complex<double> &rot_1)
{
    SymmetryKStarMember member;
    member.spatial_isym = 0;
    member.k_bz = {0.0, 0.0, 0.0};

    SymmetryKAtomRotation atom_0;
    atom_0.atom_from = 0;
    atom_0.atom_to = 1;
    atom_0.atom_type = 0;
    atom_0.lmax = 0;
    atom_0.bloch_rsh_rotations[0] = ComplexMatrix(1, 1);
    atom_0.bloch_rsh_rotations[0](0, 0) = rot_0;

    SymmetryKAtomRotation atom_1;
    atom_1.atom_from = 1;
    atom_1.atom_to = 0;
    atom_1.atom_type = 0;
    atom_1.lmax = 0;
    atom_1.bloch_rsh_rotations[0] = ComplexMatrix(1, 1);
    atom_1.bloch_rsh_rotations[0](0, 0) = rot_1;

    member.atom_rotations = {atom_0, atom_1};
    return member;
}

void test_headwing_wfc_restore_applies_atom_permutation()
{
    SymmetryContext ctx;
    ctx.set_available();
    ctx.atom_to_type = {{0, 0}, {1, 0}};
    ctx.input_coord_frac = {{0, {0.0, 0.0, 0.0}}, {1, {0.5, 0.0, 0.0}}};

    SymmetryOperation identity_operation;
    identity_operation.rotation.Identity();
    identity_operation.translation = {0.0, 0.0, 0.0};
    ctx.rspace_operations.push_back(identity_operation);

    SpeciesBasisLayout layout;
    layout.label = "X";
    layout.set({0});
    const std::vector<SpeciesBasisLayout> layouts{layout};
    const std::map<librpa_int::atom_t, size_t> atom_nw{{0, 1}, {1, 1}};
    const auto member = make_headwing_wfc_atom_swap_member(
        {2.0, 0.5}, {3.0, -0.25});

    ComplexMatrix wfc_ibz(1, 2);
    wfc_ibz(0, 0) = {0.7, -0.2};
    wfc_ibz(0, 1) = {-0.4, 0.6};

    const auto wfc_bz = librpa_int::rotate_headwing_wfc_to_kstar_member(
        ctx, member, layouts, atom_nw, {0.0, 0.0, 0.0}, wfc_ibz, nullptr);

    assert_complex_close(wfc_bz(0, 0), wfc_ibz(0, 1) * std::complex<double>{3.0, -0.25},
                         1e-12);
    assert_complex_close(wfc_bz(0, 1), wfc_ibz(0, 0) * std::complex<double>{2.0, 0.5},
                         1e-12);
}

void test_headwing_wfc_restore_applies_time_reversal()
{
    SymmetryContext ctx;
    ctx.set_available();
    ctx.atom_to_type = {{0, 0}};
    ctx.input_coord_frac = {{0, {0.0, 0.0, 0.0}}};

    SymmetryOperation identity_operation;
    identity_operation.rotation.Identity();
    identity_operation.translation = {0.0, 0.0, 0.0};
    ctx.rspace_operations.push_back(identity_operation);

    SpeciesBasisLayout layout;
    layout.label = "X";
    layout.set({0});
    const std::vector<SpeciesBasisLayout> layouts{layout};
    const std::map<librpa_int::atom_t, size_t> atom_nw{{0, 1}};

    SymmetryKStarMember member;
    member.spatial_isym = 0;
    member.k_bz = {0.0, 0.0, 0.0};
    member.time_reversal = true;
    SymmetryKAtomRotation atom;
    atom.atom_from = 0;
    atom.atom_to = 0;
    atom.atom_type = 0;
    atom.lmax = 0;
    atom.bloch_rsh_rotations[0] = ComplexMatrix(1, 1);
    atom.bloch_rsh_rotations[0](0, 0) = {0.25, 0.75};
    member.atom_rotations.push_back(atom);

    ComplexMatrix wfc_ibz(1, 1);
    wfc_ibz(0, 0) = {0.6, -0.35};

    const auto wfc_bz = librpa_int::rotate_headwing_wfc_to_kstar_member(
        ctx, member, layouts, atom_nw, {0.0, 0.0, 0.0}, wfc_ibz, nullptr);

    assert_complex_close(wfc_bz(0, 0), std::conj(wfc_ibz(0, 0)) * std::complex<double>{0.25, 0.75},
                         1e-12);
}

void test_headwing_velocity_restore_uses_inverse_spatial_route()
{
    SymmetryContext ctx;
    ctx.lattice_vectors.Identity();
    SymmetryOperation operation;
    operation.rotation = Matrix3(0.0, -1.0, 0.0,
                                 1.0,  0.0, 0.0,
                                 0.0,  0.0, 1.0);
    ctx.rspace_operations = {operation};

    SymmetryKStarMember member;
    member.spatial_isym = 0;
    std::array<ComplexMatrix, 3> velocity_ibz;
    for (auto &component : velocity_ibz) component.create(1, 1);
    velocity_ibz[0](0, 0) = {1.0, 2.0};
    velocity_ibz[1](0, 0) = {3.0, -1.0};
    velocity_ibz[2](0, 0) = {-0.5, 0.25};

    const auto velocity_bz = librpa_int::rotate_headwing_velocity_to_kstar_member(
        ctx, member, velocity_ibz, 1, false);
    assert_complex_close(velocity_bz[0](0, 0), -velocity_ibz[1](0, 0), 1e-12);
    assert_complex_close(velocity_bz[1](0, 0), velocity_ibz[0](0, 0), 1e-12);
    assert_complex_close(velocity_bz[2](0, 0), velocity_ibz[2](0, 0), 1e-12);

    member.time_reversal = true;
    const auto velocity_bz_tr = librpa_int::rotate_headwing_velocity_to_kstar_member(
        ctx, member, velocity_ibz, 1, true);
    assert_complex_close(velocity_bz_tr[0](0, 0), std::conj(velocity_ibz[1](0, 0)), 1e-12);
    assert_complex_close(velocity_bz_tr[1](0, 0), -std::conj(velocity_ibz[0](0, 0)), 1e-12);
    assert_complex_close(velocity_bz_tr[2](0, 0), -std::conj(velocity_ibz[2](0, 0)), 1e-12);
}

void test_headwing_direct_full_bz_velocity_selects_kstar_member()
{
    librpa_int::velocity_matrix_t velocity_full;
    librpa_int::initialize_velocity_matrix(velocity_full, 1, 2, 1);
    velocity_full[0][0][0](0, 0) = {1.0, 0.0};
    velocity_full[0][1][0](0, 0) = {2.0, 0.0};
    velocity_full[0][1][1](0, 0) = {3.0, 0.0};
    velocity_full[0][1][2](0, 0) = {4.0, 0.0};

    const std::vector<std::vector<int>> member_source_ik{{1}};
    const auto &velocity = librpa_int::direct_full_bz_velocity_for_kstar_member(
        velocity_full, member_source_ik, 0, 0, 0);
    assert_complex_close(velocity[0](0, 0), {2.0, 0.0}, 1e-12);
    assert_complex_close(velocity[1](0, 0), {3.0, 0.0}, 1e-12);
    assert_complex_close(velocity[2](0, 0), {4.0, 0.0}, 1e-12);
}

void test_headwing_direct_full_bz_wfc_selects_same_kstar_member()
{
    MeanField wfc_full(1, 2, 1, 1);
    auto &wfc_k0 = wfc_full.get_eigenvectors()[0][0][0];
    wfc_k0.create(1, 1);
    wfc_k0(0, 0) = {1.0, 0.0};
    auto &wfc_k1 = wfc_full.get_eigenvectors()[0][0][1];
    wfc_k1.create(1, 1);
    wfc_k1(0, 0) = {0.0, 1.0};

    const std::vector<std::vector<int>> member_source_ik{{1}};
    const auto &wfc = librpa_int::direct_full_bz_wfc_for_kstar_member(
        wfc_full, member_source_ik, 0, 0, 0, 0);
    assert_complex_close(wfc(0, 0), {0.0, 1.0}, 1e-12);
}

void test_headwing_direct_full_bz_wfc_local_block(const BlacsCtxtHandler &blacs_h)
{
    ComplexMatrix wfc_full(5, 7);
    for (int iband = 0; iband != wfc_full.nr; ++iband)
        for (int iao = 0; iao != wfc_full.nc; ++iao)
            wfc_full(iband, iao) = {100.0 * iband + iao, iband - 0.1 * iao};

    ArrayDesc desc_wfc(blacs_h);
    desc_wfc.init(7, 5, 2, 3, 0, 0);
    const auto wfc_local = librpa_int::localize_direct_full_bz_wfc(wfc_full, desc_wfc);
    assert(wfc_local.nr == desc_wfc.n_loc());
    assert(wfc_local.nc == desc_wfc.m_loc());
    for (int jloc = 0; jloc != desc_wfc.n_loc(); ++jloc)
    {
        const int iband = desc_wfc.indx_l2g_c(jloc);
        for (int iloc = 0; iloc != desc_wfc.m_loc(); ++iloc)
        {
            const int iao = desc_wfc.indx_l2g_r(iloc);
            assert_complex_close(wfc_local(jloc, iloc), wfc_full(iband, iao), 1e-12);
        }
    }
}

RI::Tensor<double> make_single_value_tensor(const double value)
{
    auto data = std::make_shared<std::valarray<double>>(1);
    (*data)[0] = value;
    return RI::Tensor<double>({1UL, 1UL, 1UL}, data);
}

void compare_local_blacs_matrices(
    const std::pair<ArrayDesc, matrix_m<std::complex<double>>> &actual,
    const std::pair<ArrayDesc, matrix_m<std::complex<double>>> &expected,
    const double tolerance)
{
    assert(actual.first.m() == expected.first.m());
    assert(actual.first.n() == expected.first.n());
    assert(actual.first.m_loc() == expected.first.m_loc());
    assert(actual.first.n_loc() == expected.first.n_loc());
    for (int i = 0; i != actual.first.m_loc(); ++i)
    {
        for (int j = 0; j != actual.first.n_loc(); ++j)
        {
            assert_complex_close(actual.second(i, j), expected.second(i, j), tolerance);
        }
    }
}

void test_kblacs_transform_with_restored_wfc_matches_full_bz_atom_permutation(
    const BlacsCtxtHandler &blacs_h)
{
    SymmetryContext ctx;
    ctx.set_available();
    ctx.atom_to_type = {{0, 0}, {1, 0}};
    ctx.input_coord_frac = {{0, {0.0, 0.0, 0.0}}, {1, {0.5, 0.0, 0.0}}};

    SymmetryOperation identity_operation;
    identity_operation.rotation.Identity();
    identity_operation.translation = {0.0, 0.0, 0.0};
    ctx.rspace_operations.push_back(identity_operation);

    SpeciesBasisLayout layout;
    layout.label = "X";
    layout.set({0});
    const std::vector<SpeciesBasisLayout> layouts{layout};
    const std::map<librpa_int::atom_t, size_t> atom_nw{{0, 1}, {1, 1}};

    auto member = make_headwing_wfc_atom_swap_member({0.0, 1.0}, {-1.0, 0.0});
    member.k_bz = {0.5, 0.0, 0.0};

    MeanField mf_ibz(1, 1, 2, 2, 1);
    auto &wfc_ibz = mf_ibz.get_eigenvectors()[0][0][0];
    wfc_ibz.create(2, 2);
    wfc_ibz(0, 0) = {0.7, -0.2};
    wfc_ibz(0, 1) = {-0.4, 0.6};
    wfc_ibz(1, 0) = {0.3, 0.5};
    wfc_ibz(1, 1) = {-0.8, -0.1};

    const auto wfc_bz = librpa_int::rotate_headwing_wfc_to_kstar_member(
        ctx, member, layouts, atom_nw, {0.0, 0.0, 0.0}, wfc_ibz, &member.k_bz);

    MeanField mf_full(1, 1, 2, 2, 1);
    auto &wfc_full = mf_full.get_eigenvectors()[0][0][0];
    wfc_full.create(2, 2);
    for (int ib = 0; ib != 2; ++ib)
    {
        for (int iao = 0; iao != 2; ++iao)
        {
            wfc_full(ib, iao) = wfc_bz(ib, iao);
        }
    }

    librpa_int::velocity_matrix_t velocity;
    librpa_int::initialize_velocity_matrix(velocity, 1, 1, 2);
    AtomicBasis basis_wfc(std::vector<size_t>{1, 1});
    AtomicBasis basis_abf(std::vector<size_t>{1, 1});
    PeriodicBoundaryData pbc;
    const std::vector<Vector3_Order<double>> kfrac_ibz{{0.0, 0.0, 0.0}};
    const std::vector<double> omega{0.5};

    diele_func df_ibz(mf_ibz, velocity, kfrac_ibz, basis_wfc, basis_abf, omega, 2, 2, 1, 1,
                      pbc, librpa_int::global::mpi_comm_global_h, blacs_h);
    diele_func df_full(mf_full, velocity, {member.k_bz}, basis_wfc, basis_abf, omega, 2, 2, 1, 1,
                       pbc, librpa_int::global::mpi_comm_global_h, blacs_h);

    std::map<int, std::map<librpa_int::libri_types<int, int>::TAC, RI::Tensor<double>>> Cs_IJ;
    Cs_IJ[0][{0, {0, 0, 0}}] = make_single_value_tensor(1.0);
    Cs_IJ[0][{1, {0, 0, 0}}] = make_single_value_tensor(-0.35);
    Cs_IJ[0][{1, {1, 0, 0}}] = make_single_value_tensor(0.42);

    std::vector<std::vector<const ComplexMatrix *>> restored_wfc_ptrs(
        1, std::vector<const ComplexMatrix *>(1, &wfc_bz));
    const auto restored = df_ibz.transform_Cs2mnk_kblacs(
        0, 0, Cs_IJ, blacs_h, member.k_bz, &restored_wfc_ptrs);
    const auto full = df_full.transform_Cs2mnk_kblacs(
        0, 0, Cs_IJ, blacs_h, member.k_bz);

    compare_local_blacs_matrices(restored, full, 1e-12);
}

void test_kblacs_transform_matches_original_transform(const BlacsCtxtHandler &blacs_h)
{
    const int nprocs = librpa_int::get_mpi_size(MPI_COMM_WORLD);
    KPointBlacsProcessShape shape(1, nprocs, true);
    KPointBlacsParallelContext kctx(shape, MPI_COMM_WORLD, 1);
    const auto desc_wfc = kctx.create_array_desc(2, 2, 2, 2);

    MeanField mf(1, 1, 2, 2, 1);
    auto &wfc = mf.get_eigenvectors()[0][0][0];
    wfc.create(2, 2);
    wfc(0, 0) = {0.7, 0.2};
    wfc(0, 1) = {-0.3, 0.4};
    wfc(1, 0) = {0.5, -0.1};
    wfc(1, 1) = {0.9, 0.3};

    librpa_int::velocity_matrix_t velocity;
    librpa_int::initialize_velocity_matrix(velocity, 1, 1, 2);
    AtomicBasis basis_wfc({2});
    AtomicBasis basis_abf({1});
    PeriodicBoundaryData pbc;
    const std::vector<Vector3_Order<double>> kfrac{{0.25, 0.0, 0.0}};
    const std::vector<double> omega{0.5};

    diele_func df(mf, velocity, kfrac, basis_wfc, basis_abf, omega, 2, 2, 1, 1, pbc,
                  librpa_int::global::mpi_comm_global_h, blacs_h, true, &kctx, &desc_wfc);

    auto tensor_R0 = std::make_shared<std::valarray<double>>(4);
    (*tensor_R0)[0] = 1.0;
    (*tensor_R0)[1] = 0.2;
    (*tensor_R0)[2] = -0.4;
    (*tensor_R0)[3] = 0.8;
    auto tensor_R1 = std::make_shared<std::valarray<double>>(4);
    (*tensor_R1)[0] = 0.3;
    (*tensor_R1)[1] = -0.1;
    (*tensor_R1)[2] = 0.2;
    (*tensor_R1)[3] = 0.5;
    std::map<int, std::map<librpa_int::libri_types<int, int>::TAC, RI::Tensor<double>>> Cs_IJ;
    Cs_IJ[0][{0, {0, 0, 0}}] = RI::Tensor<double>({1UL, 2UL, 2UL}, tensor_R0);
    Cs_IJ[0][{0, {1, 0, 0}}] = RI::Tensor<double>({1UL, 2UL, 2UL}, tensor_R1);

    librpa_int::Cs_LRI Cs_data;
    Cs_data.use_libri = true;
    Cs_data.data_libri = Cs_IJ;
    const auto targets = librpa_int::build_headwing_full_bz_fourier_targets(kfrac);
    const auto Cs_IJ_k =
        librpa_int::fourier_headwing_cs_to_ijk(Cs_data, basis_wfc, basis_abf, targets);

    auto original = df.transform_Cs2mnk(0, 0, Cs_IJ);
    auto kblacs = df.transform_Cs2mnk_kblacs(0, 0, Cs_IJ, kctx.blacs_h, kfrac[0]);
    auto prefourier = df.transform_Cs2mnk_kblacs(0, 0, 0, Cs_IJ_k, kctx.blacs_h);

    assert(original.first.m() == kblacs.first.m());
    assert(original.first.n() == kblacs.first.n());
    assert(original.first.m_loc() == kblacs.first.m_loc());
    assert(original.first.n_loc() == kblacs.first.n_loc());
    for (int i = 0; i != original.first.m_loc(); ++i)
    {
        for (int j = 0; j != original.first.n_loc(); ++j)
        {
            assert_complex_close(kblacs.second(i, j), original.second(i, j), 1e-12);
            assert_complex_close(prefourier.second(i, j), kblacs.second(i, j), 1e-12);
        }
    }

    std::vector<std::vector<ComplexMatrix>> override_storage(1);
    override_storage[0].resize(1);
    std::vector<std::vector<const ComplexMatrix *>> override_ptrs(1);
    override_ptrs[0].assign(1, nullptr);
    if (desc_wfc.is_src())
    {
        auto &override_wfc = override_storage[0][0];
        override_wfc.create(2, 2);
        override_wfc(0, 0) = {0.2, -0.4};
        override_wfc(0, 1) = {0.6, 0.1};
        override_wfc(1, 0) = {-0.5, 0.3};
        override_wfc(1, 1) = {0.4, -0.2};
        override_ptrs[0][0] = &override_wfc;
    }
    const auto real_override = df.transform_Cs2mnk_kblacs(
        0, 0, Cs_IJ, kctx.blacs_h, kfrac[0], &override_ptrs);
    const auto prefourier_override = df.transform_Cs2mnk_kblacs(
        0, 0, 0, Cs_IJ_k, kctx.blacs_h, &override_ptrs);
    int override_changed_local = 0;
    for (int i = 0; i != real_override.first.m_loc(); ++i)
    {
        for (int j = 0; j != real_override.first.n_loc(); ++j)
        {
            assert_complex_close(prefourier_override.second(i, j), real_override.second(i, j),
                                 1e-12);
            if (std::abs(real_override.second(i, j) - kblacs.second(i, j)) > 1e-12)
                override_changed_local = 1;
        }
    }
    int override_changed = 0;
    MPI_Allreduce(&override_changed_local, &override_changed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    assert(override_changed == 1);
}

void test_kblacs_transform_uses_rectangular_opt_128_blocks()
{
    constexpr int n_basis = 320;
    constexpr int n_states = 300;
    constexpr int n_ao_Mu = 160;

    const int nprocs = librpa_int::get_mpi_size(MPI_COMM_WORLD);
    KPointBlacsProcessShape shape(1, nprocs, true);
    KPointBlacsParallelContext kctx(shape, MPI_COMM_WORLD, 1);
    const auto desc_wfc_full =
        kctx.create_array_desc(n_basis, n_states, n_basis, n_states);
    const int block_ao =
        librpa_int::get_capped_blacs_block_size(
            n_basis, librpa_int::wfc_gemm_block_size_opt, kctx.blacs_h);
    const int block_band =
        librpa_int::get_capped_blacs_block_size(
            n_states, librpa_int::wfc_gemm_block_size_opt, kctx.blacs_h);
    const auto desc_wfc =
        kctx.create_array_desc(n_basis, n_states, block_ao, block_band);

    MeanField mf(1, 1, n_states, n_basis, 1);
    if (desc_wfc_full.is_src())
    {
        auto &wfc = mf.get_eigenvectors()[0][0][0];
        wfc.create(n_states, n_basis);
        wfc.zero_out();
        for (int ib = 0; ib != n_states; ++ib) wfc(ib, ib) = {1.0, 0.0};
    }
    librpa_int::redistribute_meanfield_eigvecs_kblacs(
        mf, kctx, desc_wfc_full, desc_wfc, "test headwing");

    librpa_int::velocity_matrix_t velocity;
    librpa_int::initialize_velocity_matrix(velocity, 1, 1, n_states);
    AtomicBasis basis_wfc(std::vector<std::size_t>{n_ao_Mu, n_basis - n_ao_Mu});
    AtomicBasis basis_abf(std::vector<std::size_t>{1, 1});
    PeriodicBoundaryData pbc;
    const std::vector<Vector3_Order<double>> kfrac{{0.0, 0.0, 0.0}};
    const std::vector<double> omega{0.5};

    diele_func df(mf, velocity, kfrac, basis_wfc, basis_abf, omega, n_basis, n_states, 1,
                  2, pbc, librpa_int::global::mpi_comm_global_h, kctx.blacs_h, true, &kctx,
                  &desc_wfc);

    auto tensor_data = std::make_shared<std::valarray<std::complex<double>>>(
        std::complex<double>{0.0, 0.0}, static_cast<std::size_t>(n_ao_Mu) * n_ao_Mu);
    for (int i = 0; i != n_ao_Mu; ++i)
        (*tensor_data)[static_cast<std::size_t>(i) * n_ao_Mu + i] = {1.0, 0.0};

    librpa_int::HeadwingCsIJKMap Cs_IJ_k;
    Cs_IJ_k[0][{0, 0}] = RI::Tensor<std::complex<double>>(
        {1UL, static_cast<std::size_t>(n_ao_Mu), static_cast<std::size_t>(n_ao_Mu)},
        tensor_data);

    const auto transformed =
        df.transform_Cs2mnk_kblacs(0, 0, 0, Cs_IJ_k, kctx.blacs_h);
    const auto &desc = transformed.first;
    const auto &matrix = transformed.second;
    assert(desc.m() == n_states && desc.n() == n_states);
    assert(desc.mb() == librpa_int::wfc_gemm_block_size_opt &&
           desc.nb() == librpa_int::wfc_gemm_block_size_opt);

    for (int ilo = 0; ilo != desc.m_loc(); ++ilo)
    {
        const int i = desc.indx_l2g_r(ilo);
        for (int jlo = 0; jlo != desc.n_loc(); ++jlo)
        {
            const int j = desc.indx_l2g_c(jlo);
            const std::complex<double> expected =
                i == j && i < n_ao_Mu ? std::complex<double>{2.0, 0.0}
                                       : std::complex<double>{0.0, 0.0};
            assert_complex_close(matrix(ilo, jlo), expected, 1e-12);
        }
    }
}

void test_transform_Cs2mnk_can_keep_spin_channels_separate(const BlacsCtxtHandler &blacs_h)
{
    const int nprocs = librpa_int::get_mpi_size(MPI_COMM_WORLD);
    KPointBlacsProcessShape shape(1, nprocs, true);
    KPointBlacsParallelContext kctx(shape, MPI_COMM_WORLD, 1);
    const auto desc_wfc = kctx.create_array_desc(2, 2, 2, 2);

    MeanField mf(2, 1, 2, 2, 1);
    auto &wfc_up = mf.get_eigenvectors()[0][0][0];
    auto &wfc_dn = mf.get_eigenvectors()[1][0][0];
    wfc_up.create(2, 2);
    wfc_dn.create(2, 2);
    wfc_up(0, 0) = {0.7, 0.2};
    wfc_up(0, 1) = {-0.3, 0.4};
    wfc_up(1, 0) = {0.5, -0.1};
    wfc_up(1, 1) = {0.9, 0.3};
    wfc_dn(0, 0) = {0.2, -0.6};
    wfc_dn(0, 1) = {0.8, 0.1};
    wfc_dn(1, 0) = {-0.4, 0.5};
    wfc_dn(1, 1) = {0.6, -0.2};

    librpa_int::velocity_matrix_t velocity;
    librpa_int::initialize_velocity_matrix(velocity, 2, 1, 2);
    AtomicBasis basis_wfc({2});
    AtomicBasis basis_abf({1});
    PeriodicBoundaryData pbc;
    const std::vector<Vector3_Order<double>> kfrac{{0.0, 0.0, 0.0}};
    const std::vector<double> omega{0.5};

    diele_func df(mf, velocity, kfrac, basis_wfc, basis_abf, omega, 2, 2, 2, 1, pbc,
                  librpa_int::global::mpi_comm_global_h, blacs_h, true, &kctx, &desc_wfc);

    auto tensor_data = std::make_shared<std::valarray<double>>(4);
    (*tensor_data)[0] = 1.0;
    (*tensor_data)[1] = 0.2;
    (*tensor_data)[2] = -0.4;
    (*tensor_data)[3] = 0.8;
    std::map<int, std::map<librpa_int::libri_types<int, int>::TAC, RI::Tensor<double>>> Cs_IJ;
    Cs_IJ[0][{0, {0, 0, 0}}] = RI::Tensor<double>({1UL, 2UL, 2UL}, tensor_data);
    librpa_int::Cs_LRI Cs_data;
    Cs_data.use_libri = true;
    Cs_data.data_libri = Cs_IJ;
    const auto Cs_IJ_k = librpa_int::fourier_headwing_cs_to_ijk(
        Cs_data, basis_wfc, basis_abf,
        librpa_int::build_headwing_full_bz_fourier_targets(kfrac));

    const auto all_spin = df.transform_Cs2mnk(0, 0, Cs_IJ);
    const auto spin_up = df.transform_Cs2mnk(0, 0, Cs_IJ, 0);
    const auto spin_dn = df.transform_Cs2mnk(0, 0, Cs_IJ, 1);
    const auto all_spin_k = df.transform_Cs2mnk_kblacs(0, 0, 0, Cs_IJ_k, kctx.blacs_h);
    const auto spin_up_k = df.transform_Cs2mnk_kblacs(0, 0, 0, Cs_IJ_k, kctx.blacs_h, nullptr, 0);
    const auto spin_dn_k = df.transform_Cs2mnk_kblacs(0, 0, 0, Cs_IJ_k, kctx.blacs_h, nullptr, 1);

    int spin_channels_differ_local = 0;
    for (int i = 0; i != all_spin.first.m_loc(); ++i)
    {
        for (int j = 0; j != all_spin.first.n_loc(); ++j)
        {
            assert_complex_close(all_spin.second(i, j), spin_up.second(i, j) + spin_dn.second(i, j),
                                 1e-12);
            assert_complex_close(all_spin_k.second(i, j),
                                 spin_up_k.second(i, j) + spin_dn_k.second(i, j), 1e-12);
            if (std::abs(spin_up_k.second(i, j) - spin_dn_k.second(i, j)) > 1e-12)
                spin_channels_differ_local = 1;
        }
    }
    int spin_channels_differ = 0;
    MPI_Allreduce(&spin_channels_differ_local, &spin_channels_differ, 1, MPI_INT, MPI_MAX,
                  MPI_COMM_WORLD);
    assert(spin_channels_differ == 1);
}

void test_head_initialization_does_not_require_coulomb_diagonalization(
    const BlacsCtxtHandler &blacs_h)
{
    MeanField mf(1, 1, 2, 1);
    mf.get_eigenvals()[0](0, 0) = -0.5;
    mf.get_eigenvals()[0](0, 1) = 0.5;
    mf.get_weight()[0](0, 0) = 2.0;
    mf.get_weight()[0](0, 1) = 0.0;

    librpa_int::velocity_matrix_t velocity;
    librpa_int::initialize_velocity_matrix(velocity, 1, 1, 2);
    for (int alpha = 0; alpha != 3; ++alpha)
    {
        velocity[0][0][alpha](1, 0) = std::complex<double>{0.1 * (alpha + 1), 0.0};
        velocity[0][0][alpha](0, 1) = std::complex<double>{0.1 * (alpha + 1), 0.0};
    }

    AtomicBasis basis_wfc({1});
    AtomicBasis basis_abf({1});
    PeriodicBoundaryData pbc;
    const std::vector<Vector3_Order<double>> kfrac{{0.0, 0.0, 0.0}};
    const std::vector<double> omega{0.5};
    const atpair_k_cplx_mat_t empty_vq;

    diele_func df(mf, velocity, kfrac, basis_wfc, basis_abf, omega, 1, 2, 1, 1, pbc,
                  librpa_int::global::mpi_comm_global_h, blacs_h);

    df.init(0.0, empty_vq);
    df.cal_head();
    assert(df.get_head_vec().size() == 1);
}

// Dense old Coulomb-basis averaged inverse dielectric reference. sqrt(V) is an
// independent fixed input; U supplies the Coulomb eigenvectors (x1 = U[:,0] and
// the rotation back to the ABF basis). Returns eps_inv in the ABF basis.
Matz coulomb_basis_eps_inv_reference(
    const Matz &U, const Matz &sqrtV, const Matz &chi0,
    const Matz &wing_mu, const Matz &head,
    const std::vector<std::array<double, 3>> &q_pts, const std::vector<double> &q_rho,
    int n_nonsingular = -1)
{
    const int n = U.nr();
    if (n_nonsingular < 0) n_nonsingular = n;
    assert(n_nonsingular > 0 && n_nonsingular <= n);
    const int nl = n_nonsingular - 1;

    // sqrtveig = sqrt(V) * U  (= U * diag(sqrt(lambda)) when consistent)
    const auto sqrtveig = sqrtV * U;

    // E_coul = I - sqrtveig^H * chi0 * sqrtveig = U^H * E_abf * U
    auto E_coul = sqrtveig.get_transpose(true) * chi0 * sqrtveig;
    E_coul *= -1.0;
    for (int i = 0; i < n; ++i) E_coul(i, i) += 1.0;

    if (nl == 0)
    {
        // No body channels: L = H and the averaged inverse is a0 * P with
        // P = x1*x1^H. Avoids a zero-dimensional LAPACK inversion.
        const auto L00 = head(0, 0), L01 = head(0, 1), L02 = head(0, 2);
        const auto L10 = head(1, 0), L11 = head(1, 1), L12 = head(1, 2);
        const auto L20 = head(2, 0), L21 = head(2, 1), L22 = head(2, 2);
        std::complex<double> a0 = 0.0;
        for (std::size_t ileb = 0; ileb < q_pts.size(); ++ileb)
        {
            const double qx = q_pts[ileb][0], qy = q_pts[ileb][1], qz = q_pts[ileb][2];
            const auto qLq = qx * (qx * L00 + qy * L01 + qz * L02) +
                             qy * (qx * L10 + qy * L11 + qz * L12) +
                             qz * (qx * L20 + qy * L21 + qz * L22);
            a0 += q_rho[ileb] / qLq;
        }
        Matz x1(n, 1, MAJOR::COL);
        for (int i = 0; i < n; ++i) x1(i, 0) = U(i, 0);
        return a0 * (x1 * x1.get_transpose(true));
    }

    // body = E_coul[1:,1:], invert it via LU
    Matz body(nl, nl, MAJOR::COL);
    for (int i = 0; i < nl; ++i)
        for (int j = 0; j < nl; ++j)
            body(i, j) = E_coul(i + 1, j + 1);
    auto body_inv = body.copy();
    std::vector<int> ipiv(static_cast<std::size_t>(nl));
    std::vector<std::complex<double>> work(static_cast<std::size_t>(nl * nl));
    int info = 0;
    librpa_int::LapackConnector::getrf_f(
        nl, nl, body_inv.ptr(), nl, ipiv.data(), info);
    assert(info == 0);
    const int lwork = nl * nl;
    librpa_int::LapackConnector::getri_f(
        nl, body_inv.ptr(), nl, ipiv.data(), work.data(), lwork, info);
    assert(info == 0);

    // wing = sqrtveig[:,1:]^H * wing_mu  (nl x 3)
    Matz wing(nl, 3, MAJOR::COL);
    for (int i = 0; i < nl; ++i)
        for (int j = 0; j < 3; ++j)
        {
            std::complex<double> s = 0.0;
            for (int k = 0; k < n; ++k)
                s += std::conj(sqrtveig(k, i + 1)) * wing_mu(k, j);
            wing(i, j) = s;
        }

    auto Lind = head - wing.get_transpose(true) * body_inv * wing;
    auto bw = body_inv * wing;
    auto wb = wing.get_transpose(true) * body_inv;

    const auto L00 = Lind(0, 0), L01 = Lind(0, 1), L02 = Lind(0, 2);
    const auto L10 = Lind(1, 0), L11 = Lind(1, 1), L12 = Lind(1, 2);
    const auto L20 = Lind(2, 0), L21 = Lind(2, 1), L22 = Lind(2, 2);

    const int nleb = static_cast<int>(q_pts.size());
    Matz eps_inv_coul(n, n, MAJOR::COL);
    eps_inv_coul.zero_out();
    for (int ileb = 0; ileb < nleb; ++ileb)
    {
        const double qx = q_pts[ileb][0], qy = q_pts[ileb][1], qz = q_pts[ileb][2];
        const auto qLq = qx * (qx * L00 + qy * L01 + qz * L02) +
                         qy * (qx * L10 + qy * L11 + qz * L12) +
                         qz * (qx * L20 + qy * L21 + qz * L22);
        const auto w = q_rho[ileb] / qLq;
        eps_inv_coul(0, 0) += w;
        for (int i = 1; i < n_nonsingular; ++i)
            for (int j = 1; j < n_nonsingular; ++j)
            {
                const auto bwq = bw(i - 1, 0) * qx + bw(i - 1, 1) * qy + bw(i - 1, 2) * qz;
                const auto qwb = qx * wb(0, j - 1) + qy * wb(1, j - 1) + qz * wb(2, j - 1);
                eps_inv_coul(i, j) += w * bwq * qwb;
            }
    }
    for (int i = 1; i < n_nonsingular; ++i)
        for (int j = 1; j < n_nonsingular; ++j)
            eps_inv_coul(i, j) += body_inv(i - 1, j - 1);

    return U * eps_inv_coul * U.get_transpose(true);
}

// Run the production ABF-space helper on distributed matrices and compare to
// the dense old Coulomb-basis reference.
void run_abf_case(const BlacsCtxtHandler &blacs_h, int n, int block_size,
                  const Matz &U, const Matz &sqrtV,
                  const Matz &chi0, const Matz &wing_mu,
                  const Matz &head,
                  const std::vector<std::array<double, 3>> &q_pts,
                  const std::vector<double> &q_rho, bool use_cholesky, double tol,
                  int n_nonsingular = -1)
{
    if (n_nonsingular < 0) n_nonsingular = n;
    const auto ref = coulomb_basis_eps_inv_reference(U, sqrtV, chi0, wing_mu, head,
                                                     q_pts, q_rho, n_nonsingular);
    // E = I - sqrt(V) * chi0 * sqrt(V)
    auto E = sqrtV * chi0 * sqrtV;
    E *= -1.0;
    for (int i = 0; i < n; ++i) E(i, i) += 1.0;

    ArrayDesc desc(blacs_h);
    desc.init(n, n, block_size, block_size, 0, 0);
    auto E_dist = init_local_mat<std::complex<double>>(desc, MAJOR::COL);
    auto sqrtV_dist = init_local_mat<std::complex<double>>(desc, MAJOR::COL);
    auto U_dist = init_local_mat<std::complex<double>>(desc, MAJOR::COL);
    fill_distributed_matrix(E_dist, desc, E);
    fill_distributed_matrix(sqrtV_dist, desc, sqrtV);
    fill_distributed_matrix(U_dist, desc, U);

    std::vector<double> qx(q_pts.size()), qy(q_pts.size()), qz(q_pts.size());
    for (std::size_t i = 0; i < q_pts.size(); ++i)
    {
        qx[i] = q_pts[i][0];
        qy[i] = q_pts[i][1];
        qz[i] = q_pts[i][2];
    }

    librpa_int::rewrite_eps_abf_space(E_dist, sqrtV_dist, U_dist, head, wing_mu, qx, qy,
                                      qz, q_rho, desc, blacs_h,
                                      static_cast<std::size_t>(n_nonsingular), 0.0,
                                      use_cholesky, false);

    double local_max_err = 0.0;
    for (int ilo = 0; ilo != desc.m_loc(); ++ilo)
    {
        const int ig = desc.indx_l2g_r(ilo);
        for (int jlo = 0; jlo != desc.n_loc(); ++jlo)
        {
            const int jg = desc.indx_l2g_c(jlo);
            local_max_err =
                std::max(local_max_err, std::abs(E_dist(ilo, jlo) - ref(ig, jg)));
        }
    }
    double max_err = 0.0;
    MPI_Allreduce(&local_max_err, &max_err, 1, MPI_DOUBLE, MPI_MAX, desc.comm());
    require_double_close(max_err, 0.0, tol);
}

void test_abf_space_wing_rewrite_matches_coulomb_basis(const BlacsCtxtHandler &blacs_h)
{
    // 3D quadrature: 6 points on the unit sphere.
    const std::vector<std::array<double, 3>> q3d{
        {1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}};
    std::vector<double> rho3d(6);
    for (auto &r : rho3d) r = (4.0 * M_PI / 6.0) / 3.0;

    // 2D quadrature: 4 points on the unit circle.
    const std::vector<std::array<double, 3>> q2d{
        {1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}};
    std::vector<double> rho2d(4);
    for (auto &r : rho2d) r = (2.0 * M_PI / 4.0) / 2.0;

    constexpr int n = 4;
    const double half = 0.5;
    const Matz U(
        {{{half, 0.0}, {half, 0.0}, {half, 0.0}, {half, 0.0}},
         {{half, 0.0}, {-half, 0.0}, {half, 0.0}, {-half, 0.0}},
         {{half, 0.0}, {half, 0.0}, {-half, 0.0}, {-half, 0.0}},
         {{half, 0.0}, {-half, 0.0}, {-half, 0.0}, {half, 0.0}}},
        MAJOR::COL);
    const std::array<double, n> lambda{{9.0, 4.0, 1.0, 0.25}};
    Matz sqrtveig(n, n, MAJOR::COL);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            sqrtveig(i, j) = U(i, j) * std::sqrt(lambda[static_cast<size_t>(j)]);
    const auto sqrtV = sqrtveig * U.get_transpose(true);

    // Negative-semidefinite chi0 = -v v^H so E (and hence M) is Hermitian
    // positive definite and Cholesky is valid.
    Matz chi0(n, n, MAJOR::COL);
    {
        const std::array<std::complex<double>, n> v{{
            {0.20, 0.0}, {0.10, 0.05}, {-0.15, 0.0}, {0.05, -0.10}}};
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                chi0(i, j) = -v[static_cast<size_t>(i)] * std::conj(v[static_cast<size_t>(j)]);
    }

    const Matz wing_mu(
        {{{0.11, 0.01}, {0.12, 0.02}, {0.13, 0.03}},
         {{0.21, 0.04}, {0.22, 0.05}, {0.23, 0.06}},
         {{0.31, 0.07}, {0.32, 0.08}, {0.33, 0.09}},
         {{0.41, 0.10}, {0.42, 0.11}, {0.43, 0.12}}},
        MAJOR::COL);
    Matz zero_wing_mu(n, 3, MAJOR::COL);
    zero_wing_mu.zero_out();

    const Matz head(
        {{{2.0, 0.0}, {0.2, 0.1}, {0.3, -0.1}},
         {{0.2, -0.1}, {2.2, 0.0}, {0.4, 0.2}},
         {{0.3, 0.1}, {0.4, -0.2}, {2.4, 0.0}}},
        MAJOR::COL);

    // LU and Cholesky, 3D and 2D, nonzero and zero wing (square grid).
    run_abf_case(blacs_h, n, 1, U, sqrtV, chi0, wing_mu, head, q3d, rho3d, false, 1e-10);
    run_abf_case(blacs_h, n, 2, U, sqrtV, chi0, wing_mu, head, q3d, rho3d, true, 1e-10);
    run_abf_case(blacs_h, n, 1, U, sqrtV, chi0, wing_mu, head, q2d, rho2d, false, 1e-10);
    run_abf_case(blacs_h, n, 2, U, sqrtV, chi0, wing_mu, head, q2d, rho2d, true, 1e-10);
    run_abf_case(blacs_h, n, 1, U, sqrtV, chi0, zero_wing_mu, head, q3d, rho3d, false, 1e-10);
    run_abf_case(blacs_h, n, 2, U, sqrtV, chi0, zero_wing_mu, head, q3d, rho3d, true, 1e-10);

    // Wing with a deliberately large component parallel to x1 in each
    // Cartesian column. The Direct-Z path (A = D*Z, no explicit projection)
    // must still match the Coulomb-basis reference because D*x1 = 0.
    {
        auto wing_mu_large_x1 = wing_mu.copy();
        for (int alpha = 0; alpha != 3; ++alpha)
            for (int mu = 0; mu != n; ++mu)
                wing_mu_large_x1(mu, alpha) += 100.0 * U(mu, 0);
        run_abf_case(blacs_h, n, 1, U, sqrtV, chi0, wing_mu_large_x1, head, q3d, rho3d,
                     false, 1e-10);
        run_abf_case(blacs_h, n, 2, U, sqrtV, chi0, wing_mu_large_x1, head, q3d, rho3d,
                     true, 1e-10);
    }

    // Rectangular (horizontal) CPU BLACS grid; exercises empty thin local
    // blocks on ranks that own no rows/columns of the thin descriptors.
    {
        BlacsCtxtHandler horizontal_blacs_h(MPI_COMM_WORLD);
        horizontal_blacs_h.init();
        horizontal_blacs_h.set_horizontal_grid();
        run_abf_case(horizontal_blacs_h, n, 1, U, sqrtV, chi0, wing_mu, head, q3d, rho3d,
                     false, 1e-10);
        run_abf_case(horizontal_blacs_h, n, 1, U, sqrtV, chi0, wing_mu, head, q3d, rho3d,
                     true, 1e-10);
    }

    // n_abf == 1 exercises the n-by-1 thin descriptors at the minimal size.
    {
        Matz U1(1, 1, MAJOR::COL);
        U1(0, 0) = {1.0, 0.0};
        Matz sqrtV1(1, 1, MAJOR::COL);
        sqrtV1(0, 0) = {3.0, 0.0};
        Matz chi01(1, 1, MAJOR::COL);
        chi01(0, 0) = {-0.04, 0.0};
        Matz wing_mu1(1, 3, MAJOR::COL);
        wing_mu1(0, 0) = {0.10, 0.01};
        wing_mu1(0, 1) = {0.11, 0.02};
        wing_mu1(0, 2) = {0.12, 0.03};
        run_abf_case(blacs_h, 1, 1, U1, sqrtV1, chi01, wing_mu1, head, q3d, rho3d, false,
                     1e-10);
    }

    // Filtered Coulomb basis: the last eigenchannel has zero eigenvalue and is
    // excluded. Both inversion paths must agree with the old reduced Coulomb-
    // basis algorithm and produce no component in the filtered subspace.
    {
        constexpr int nr = n - 1;
        Matz sqrtveig_filtered(n, n, MAJOR::COL);
        sqrtveig_filtered.zero_out();
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < nr; ++j)
                sqrtveig_filtered(i, j) =
                    U(i, j) * std::sqrt(lambda[static_cast<std::size_t>(j)]);
        const auto sqrtV_filtered = sqrtveig_filtered * U.get_transpose(true);

        run_abf_case(blacs_h, n, 1, U, sqrtV_filtered, chi0, wing_mu, head,
                     q3d, rho3d, false, 1e-10, nr);
        run_abf_case(blacs_h, n, 2, U, sqrtV_filtered, chi0, wing_mu, head,
                     q3d, rho3d, true, 1e-10, nr);
        run_abf_case(blacs_h, n, 1, U, sqrtV_filtered, chi0, zero_wing_mu, head,
                     q2d, rho2d, false, 1e-10, nr);
    }

    // Invariance: with E and sqrt(V) held fixed, changing Coulomb eigenvector
    // columns other than x1 must not change the ABF-space result.
    {
        auto U_mod = U.copy();
        for (int col = 1; col < n; ++col)
            for (int row = 0; row < n; ++row)
                U_mod(row, col) *= std::complex<double>{0.0, 2.0};

        auto E = sqrtV * chi0 * sqrtV;
        E *= -1.0;
        for (int i = 0; i < n; ++i) E(i, i) += 1.0;

        ArrayDesc desc(blacs_h);
        desc.init(n, n, 1, 1, 0, 0);
        auto E1 = init_local_mat<std::complex<double>>(desc, MAJOR::COL);
        auto E2 = init_local_mat<std::complex<double>>(desc, MAJOR::COL);
        auto sqrtV_dist = init_local_mat<std::complex<double>>(desc, MAJOR::COL);
        auto U_orig_dist = init_local_mat<std::complex<double>>(desc, MAJOR::COL);
        auto U_mod_dist = init_local_mat<std::complex<double>>(desc, MAJOR::COL);
        fill_distributed_matrix(E1, desc, E);
        fill_distributed_matrix(E2, desc, E);
        fill_distributed_matrix(sqrtV_dist, desc, sqrtV);
        fill_distributed_matrix(U_orig_dist, desc, U);
        fill_distributed_matrix(U_mod_dist, desc, U_mod);

        std::vector<double> qx(6), qy(6), qz(6);
        for (std::size_t i = 0; i < q3d.size(); ++i)
        {
            qx[i] = q3d[i][0];
            qy[i] = q3d[i][1];
            qz[i] = q3d[i][2];
        }

        librpa_int::rewrite_eps_abf_space(E1, sqrtV_dist, U_orig_dist, head, wing_mu, qx, qy,
                                          qz, rho3d, desc, blacs_h,
                                          static_cast<std::size_t>(n), 0.0, false, false);
        librpa_int::rewrite_eps_abf_space(E2, sqrtV_dist, U_mod_dist, head, wing_mu, qx, qy,
                                          qz, rho3d, desc, blacs_h,
                                          static_cast<std::size_t>(n), 0.0, false, false);

        double local_max_err = 0.0;
        for (int ilo = 0; ilo != desc.m_loc(); ++ilo)
            for (int jlo = 0; jlo != desc.n_loc(); ++jlo)
                local_max_err =
                    std::max(local_max_err, std::abs(E1(ilo, jlo) - E2(ilo, jlo)));
        double max_err = 0.0;
        MPI_Allreduce(&local_max_err, &max_err, 1, MPI_DOUBLE, MPI_MAX, desc.comm());
        require_double_close(max_err, 0.0, 1e-10);
    }
}
}  // namespace

int main(int argc, char *argv[])
{
    int provided = 0;
    MPI_Init_thread(&argc, &argv, LIBRPA_MPI_THREAD_LEVEL, &provided);
    librpa_int::global::init_global_mpi(MPI_COMM_WORLD);
    librpa_int::global::init_global_io(false, "stdout", false);

    {
        BlacsCtxtHandler blacs_h(MPI_COMM_WORLD);
        blacs_h.init();
        blacs_h.set_square_grid();

        test_headwing_body_inverse_uses_identity_solve(blacs_h);
        test_replace_rpa_response_headwing_replaces_only_singular_channels(blacs_h);
        test_gamma_head_rank_one_matches_coulomb_basis_overwrite(blacs_h);
        test_gamma_head_rank_one_handles_empty_local_blocks(blacs_h);
        test_rspace_symmetry_requires_complete_band_space();
        test_kpoint_coordinate_mapping_selects_active_klist_from_full_source();
        test_kstar_velocity_mapping_preserves_member_order_and_periodic_gauge();
        test_replace_rpa_response_head_only_keeps_numeric_wings(blacs_h);
        test_rpa_trace_log_average_uses_directional_head_and_wing();
        test_rpa_headwing_regular_body_start_channel();
        test_rpa_headwing_gamma_cell_volume_uses_reciprocal_lattice();
        test_rpa_chi0v_wing_desc_uses_global_rows(blacs_h);
        test_headwing_spin_weights();
        test_wing_cartesian_gram_is_invariant_under_row_phases();
        test_velocity_matrix_initialization();
        test_headwing_local_kpoints_prefers_kpoint_blacs_context();
        test_headwing_world_fourier_uses_all_R_blocks_at_nonzero_k();
        test_headwing_symmetry_fourier_target_ids_are_deterministic();
        test_headwing_ijk_redistribution_is_owner_group_local();
        test_accumulate_wing_mu_for_pair_matches_original_formula();
        test_headwing_wfc_restore_applies_atom_permutation();
        test_headwing_wfc_restore_applies_time_reversal();
        test_headwing_velocity_restore_uses_inverse_spatial_route();
        test_headwing_direct_full_bz_velocity_selects_kstar_member();
        test_headwing_direct_full_bz_wfc_selects_same_kstar_member();
        test_headwing_direct_full_bz_wfc_local_block(blacs_h);
        test_kblacs_transform_with_restored_wfc_matches_full_bz_atom_permutation(blacs_h);
        test_kblacs_transform_matches_original_transform(blacs_h);
        test_kblacs_transform_uses_rectangular_opt_128_blocks();
        test_transform_Cs2mnk_can_keep_spin_channels_separate(blacs_h);
        test_head_initialization_does_not_require_coulomb_diagonalization(blacs_h);
        test_abf_space_wing_rewrite_matches_coulomb_basis(blacs_h);
    }

    librpa_int::global::finalize_global_io();
    librpa_int::global::finalize_global_mpi();
    MPI_Finalize();
    return 0;
}
