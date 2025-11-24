#include "task_test.h"

#include "../src/mpi/envs_mpi.h"
#include "../src/utils/envs_io.h"
#include "../src/utils/profiler.h"
#include "../src/math/matrix_m.h"

#include "read_data.h"
#include "driver_params.h"

void task_test()
{
    using librpa_int::envs::mpi_comm_global_h;
    using librpa_int::envs::ofs_myid;
    using librpa_int::utils::lib_printf;

    Profiler::start("test");

    int n_basis_band, n_states_band, n_spin_band;
    std::vector<Vector3_Order<double>> kfrac_band = read_band_kpath_info(
        driver_params.input_dir + "band_kpath_info", n_basis_band, n_states_band, n_spin_band);
    if (mpi_comm_global_h.is_root())
    {
        std::cout << "Band k-points to compute:\n";
        for (int ik = 0; ik < kfrac_band.size(); ik++)
        {
            const auto &k = kfrac_band[ik];
            lib_printf("%5d %12.7f %12.7f %12.7f\n", ik + 1, k.x, k.y, k.z);
        }
    }
    mpi_comm_global_h.barrier();

    auto meanfield_band = read_meanfield_band(driver_params.input_dir,
            n_basis_band, n_states_band, n_spin_band, kfrac_band.size());

    assert(n_basis_band == n_states_band);

    // print eigenvectors
    for (int ispin = 0; ispin < n_spin_band; ispin++)
    {
        for (int ik = 0; ik < kfrac_band.size(); ik++)
        {
            const matrix_m<complex<double>> vec(n_basis_band, n_states_band,
                    meanfield_band.get_eigenvectors()[ispin][ik].c, MAJOR::ROW);
            stringstream s;
            s << "bandv_" << ispin << "_" << ik << ".csc";
            write_matrix_elsi_csc(vec, s.str());
        }
    }

    Profiler::cease("test");
}
