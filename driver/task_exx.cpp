#include "task_exx.h"

// other driver headers
#include "read_data.h"

// src headers
#include "app_exx.h"
#include "constants.h"
#include "meanfield.h"
#include "envs_mpi.h"
#include "utils_io.h"
#include "pbc.h"

// debug headers
// #include "stl_io_helper.h"

void task_exx()
{
    using LIBRPA::envs::mpi_comm_global_h;
    using LIBRPA::utils::lib_printf;

    // Load the cut coulomb data
    read_Vq_full("./", "coulomb_cut_", true);

    // NOTE: since the internal object does not support parallel calculation with different k-points on each process
    // we compute all and only let master process print
    const int i_state_low = 0;
    const int n_kpoints = meanfield.get_n_kpoints();
    const int i_state_high = meanfield.get_n_bands();
    const int n_states_calc = i_state_high - i_state_low;

    auto exx_ks = LIBRPA::app::compute_exx_orbital_energy_(i_state_low, i_state_high, -1, nullptr);

    if (exx_ks.size() > 0 && mpi_comm_global_h.is_root())
    {
        for (int isp = 0; isp != meanfield.get_n_spins(); isp++)
        {
            lib_printf("Spin channel %1d\n", isp+1);
            for (int ik = 0; ik != meanfield.get_n_kpoints(); ik++)
            {
                cout << "k-point " << ik + 1 << ": " << kfrac_list[ik] << endl;
                lib_printf("%-4s  %-10s  %-10s\n", "Band", "e_exx (Ha)", "e_exx (eV)");
                for (int ib = i_state_low; ib != i_state_high; ib++)
                {
                    const int index = isp * n_kpoints * n_states_calc + ik * n_states_calc + ib - i_state_low;
                    const auto e = exx_ks[index];
                    lib_printf("%4d  %10.5f  %10.5f\n", ib + 1, e, HA2EV * e);
                }
                lib_printf("\n");
            }
        }
    }
}
