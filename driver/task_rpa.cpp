#include "task_rpa.h"

#include "../src/api/app_rpa.h"
#include "../src/mpi/global_mpi.h"
#include "../src/core/ri.h"

void task_rpa()
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
    std::complex<double> corr;
    std::vector<std::complex<double>> corr_irk(n_irk_points);

    librpa_int::app::get_rpa_correlation_energy_(corr, corr_irk);

    if (mpi_comm_global_h.is_root())
    {
        lib_printf("RPA correlation energy (Hartree)\n");
        lib_printf("| Weighted contribution from each k:\n");

        for (int i_irk = 0; i_irk < n_irk_points; i_irk++)
        {
            cout << "| " << irk_points[i_irk] << ": " << corr_irk[i_irk] << endl;
        }
        lib_printf("| Total EcRPA: %18.9f\n", corr.real());
        if (std::abs(corr.imag()) > 1.e-3)
            lib_printf("Warning: considerable imaginary part of EcRPA = %f\n", corr.imag());
    }
}
