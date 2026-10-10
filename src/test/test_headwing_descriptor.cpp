#include <algorithm>
#include <complex>
#include <iostream>
#include <stdexcept>

#include "../core/dielecmodel.h"
#include "../io/global_io.h"
#include "../math/utils_matrix_m_mpi.h"

using namespace librpa_int;

namespace
{
void check_layouts(const BlacsCtxtHandler &blacs)
{
    constexpr int nabf = 192;
    MeanField mf(1, 1, 2, 2, 1);
    velocity_matrix_t velocity;
    AtomicBasis basis_wfc({2});
    AtomicBasis basis_abf({nabf});
    PeriodicBoundaryData pbc;
    const std::vector<Vector3_Order<double>> kfrac{{0.0, 0.0, 0.0}};
    // No frequencies are needed to exercise the production descriptor guard.
    const std::vector<double> omega;
    diele_func df(mf, velocity, kfrac, basis_wfc, basis_abf, omega, 2, 2, 1, nabf, pbc,
                  global::mpi_comm_global_h, blacs);

    for (const bool gw_layout : {false, true})
    {
        ArrayDesc seed(blacs);
        if (gw_layout)
            seed.init_1b1p(nabf, nabf, 0, 0);
        else
            seed.init_square_blk(nabf, nabf, 0, 0);
        const int block = std::min(128, seed.nb());
        ArrayDesc producer(blacs);
        producer.init(nabf, nabf, block, block, 0, 0);
        auto matrix = init_local_mat<std::complex<double>>(producer, MAJOR::COL);
        df.wing_mu_to_lambda(matrix, producer, nabf);

        // One bad rank must make every rank reject before entering collective BLAS.
        const auto major = blacs.myid == 0 ? MAJOR::ROW : MAJOR::COL;
        auto row_major = init_local_mat<std::complex<double>>(producer, major);
        bool rejected = false;
        try
        {
            df.wing_mu_to_lambda(row_major, producer, nabf);
        }
        catch (const std::logic_error &)
        {
            rejected = true;
        }
        if (!rejected) throw std::runtime_error("row-major Coulomb matrix was accepted");

        ArrayDesc wrong_size(blacs);
        wrong_size.init(nabf - 1, nabf, block, block, 0, 0);
        auto undersized = init_local_mat<std::complex<double>>(wrong_size, MAJOR::COL);
        rejected = false;
        try
        {
            df.wing_mu_to_lambda(undersized, wrong_size, nabf);
        }
        catch (const std::logic_error &)
        {
            rejected = true;
        }
        if (!rejected) throw std::runtime_error("wrong global Coulomb dimensions were accepted");
    }
}
}  // namespace

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    global::init_global_mpi(MPI_COMM_WORLD);
    global::init_global_io(false, "stdout", false);
    int size = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    try
    {
        for (const bool transpose : {false, true})
        {
            BlacsCtxtHandler blacs(MPI_COMM_WORLD);
            blacs.init();
            blacs.set_grid(transpose ? 1 : size, transpose ? size : 1);
            check_layouts(blacs);
        }
    }
    catch (const std::exception &error)
    {
        std::cerr << error.what() << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 2);
    }
    global::finalize_global_io();
    global::finalize_global_mpi();
    MPI_Finalize();
    return 0;
}
