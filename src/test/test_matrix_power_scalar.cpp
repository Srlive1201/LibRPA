#include <cmath>
#include <iostream>
#include <stdexcept>

#include "../api/dataset.h"
#include "../math/utils_matrix_m_mpi.h"
#include "../mpi/global_mpi.h"

using namespace librpa_int;

namespace
{
void close(std::complex<double> actual, std::complex<double> expected)
{
    if (!std::isfinite(actual.real()) || !std::isfinite(actual.imag()) ||
        std::abs(actual - expected) > 2e-14)
        throw std::runtime_error("scalar distributed matrix power mismatch");
}

void check(double eigenvalue, double power, double threshold, double expected)
{
    Dataset ds{MPI_COMM_WORLD};
    // Use the raw context constructor: the operation needs only the BLACS grid.
    ArrayDesc ad_a(ds.blacs_h.ictxt), ad_z(ds.blacs_h.ictxt);
    // Exercise redistribution from a nonzero owner and back to a different Z owner.
    ad_a.init(1, 1, 1, 1, ds.blacs_h.nprows - 1, ds.blacs_h.npcols - 1);
    ad_z.init(1, 1, 1, 1, 0, 0);
    auto a = init_local_mat<std::complex<double>>(ad_a, MAJOR::COL);
    auto z = init_local_mat<std::complex<double>>(ad_z, MAJOR::COL);
    if (ad_a.is_src()) a(0, 0) = eigenvalue;
    double w = 0;
    std::size_t filtered = 0;
    const auto scaled = power_hemat_blacs(a, ad_a, z, ad_z, filtered, &w, power, threshold);
    close(w, eigenvalue);
    if (filtered != std::size_t(eigenvalue < threshold))
        throw std::runtime_error("scalar distributed eigenvalue filtering mismatch");
    if (ad_a.is_src()) close(a(0, 0), expected);
    if (ad_z.is_src())
    {
        close(z(0, 0), 1);
        close(scaled(0, 0), expected);
    }
}
}  // namespace

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    global::init_global_mpi(MPI_COMM_WORLD);
    global::init_global_io();
    int status = 0;
    try
    {
        check(2, 0.5, 0, std::sqrt(2.0));
        check(-0.4, 0.5, 0, 0);
        check(4, -1, 0, 0.25);
        check(-0.4, 2, -1, 0.16);
        if (global::mpi_comm_global_h.is_root()) std::cout << "scalar matrix powers passed\n";
    }
    catch (const std::exception &error)
    {
        std::cerr << error.what() << '\n';
        status = 1;
    }
    global::finalize_global_io();
    global::finalize_global_mpi();
    MPI_Finalize();
    return status;
}
