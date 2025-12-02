#include "librpa.hpp"

#include <iostream>

#include "../task.h"
#include "../driver.h"
#include "../../src/mpi/global_mpi.h"
#include "../../src/io/global_io.h"

// #include "../../src/io/stl_io_helper.h"

void driver::task_rpa()
{
    using namespace librpa_int;
    using librpa_int::global::mpi_comm_global_h;
    using librpa_int::global::lib_printf;

    // Using public API.
    // std::vector<double> temp_corr(2);
    // std::vector<double> temp_corr_irk(2 * n_irk_points);
    // get_rpa_correlation_energy(temp_corr.data(), temp_corr_irk.data());
    // std::complex<double> corr(temp_corr[0], temp_corr[1]);
    // auto dp = reinterpret_cast<std::complex<double>*>(temp_corr_irk.data());
    // std::vector<std::complex<double>> corr_irk(dp, dp + n_irk_points);

    // Using the internal work function, avoid copying between vector and double*
    double corr = 0.0;
    std::vector<std::complex<double>> corr_irk(n_ibz_kpoints);

    corr = h.get_rpa_correlation_energy(driver::opts, corr_irk);

    mpi_comm_global_h.barrier();
    if (mpi_comm_global_h.is_root())
    {
        lib_printf("RPA correlation energy (Hartree)\n");
        lib_printf("| Weighted contribution from each k:\n");

        for (int i_irk = 0; i_irk < n_ibz_kpoints; i_irk++)
        {
            std::cout << "| " << ibz_kpoints[i_irk] << ": " << corr_irk[i_irk] << std::endl;
        }
        lib_printf("| Total EcRPA: %18.9f\n", corr);
        for (int i_irk = 0; i_irk < n_ibz_kpoints; i_irk++)
        {
            const auto &im = corr_irk[i_irk].imag();
            if (std::abs(im) > 1.e-3)
                lib_printf("Warning: considerable imaginary part of EcRPA = %f\n at IBZ k-point %d\n", im, i_irk + 1);
        }
    }
    mpi_comm_global_h.barrier();
}
