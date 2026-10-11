#include <array>
#include <cmath>
#include <complex>
#include <iostream>
#include <vector>

#include "../core/epsilon.h"
#include "../math/utils_matrix_m_mpi.h"
#include "../utils/constants.h"

using namespace librpa_int;

static bool check_logdet(const BlacsCtxtHandler& blacs,
                         const std::array<std::complex<double>, 9>& values,
                         const std::complex<double>& expected, const char* label,
                         const int block_size)
{
    ArrayDesc desc(blacs);
    desc.init(3, 3, block_size, block_size, 0, 0);
    auto matrix = init_local_mat<std::complex<double>>(desc, MAJOR::COL);
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
        {
            const int row = desc.indx_g2l_r(i), col = desc.indx_g2l_c(j);
            if (row >= 0 && col >= 0) matrix(row, col) = values[3 * i + j];
        }
    std::vector<int> pivots(desc.m_loc() + desc.mb());
    int info = 0;
    const auto actual = compute_pi_det_blacs_2d(matrix, desc, pivots.data(), info);
    const double real_error = std::abs(actual.real() - expected.real());
    const double phase_error = std::abs(std::remainder(actual.imag() - expected.imag(), TWO_PI));
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0)
        std::cout << label << " block_size=" << block_size << " actual=" << actual
                  << " expected=" << expected << " info=" << info << std::endl;
    return info == 0 && std::isfinite(real_error) && std::isfinite(phase_error)
           && real_error < 1e-12 && phase_error < 1e-12;
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    bool passed = true;
    {
        BlacsCtxtHandler blacs(MPI_COMM_WORLD);
        blacs.init();
        blacs.set_square_grid();
        const auto z = std::polar(1.0, TWO_PI / 6.0);
        const auto t = std::polar(1.0, TWO_PI * 0.4);
        for (const int block_size : {1, 2, 3})
        {
            // One row swap and z^3=-1 give det=+1, not a phase of pi.
            passed &= check_logdet(blacs, {0., z, 0., z, 0., 0., 0., 0., z},
                                   {0., 0.}, "positive_det_with_complex_pivots", block_size);
            passed &= check_logdet(blacs, {-1., 0., 0., 0., 2., 0., 0., 0., 3.},
                                   {std::log(6.0), TWO_PI / 2.0}, "negative_determinant", block_size);
            passed &= check_logdet(blacs, {t, 0., 0., 0., t, 0., 0., 0., t},
                                   {0., TWO_PI * 0.2}, "complex_determinant", block_size);
            passed &= check_logdet(blacs, {2., 0., 0., 0., 3., 0., 0., 0., 5.},
                                   {std::log(30.0), 0.}, "positive_diagonal", block_size);
        }
    }
    MPI_Finalize();
    return passed ? 0 : 1;
}
