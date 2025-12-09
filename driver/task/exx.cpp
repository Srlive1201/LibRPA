// driver headers
#include "../task.h"
#include "../driver.h"
#include "../read_data.h"

// src headers
#include "../../src/utils/profiler.h"
#include "../../src/mpi/global_mpi.h"
#include "../../src/utils/constants.h"
#include "../../src/io/global_io.h"

#include "../../src/api/instance_manager.h"
// debug headers
// #include "stl_io_helper.h"

void driver::task_exx()
{
    using librpa_int::global::mpi_comm_global_h;
    using librpa_int::global::lib_printf;

    using std::cout;
    using std::endl;
    using std::map;
    using namespace librpa_int::global;
    using librpa_int::global::lib_printf;
    using librpa_int::cplxdb;
    using librpa_int::HA2EV;

    profiler.start("exx", "Exact-exchange only calculation");

    profiler.start("read_vq_cut", "Load truncated Coulomb");

    // Prepare input specific to G0W0
    auto routing = opts.parallel_routing;
    if (routing == LibrpaParallelRouting::AUTO) routing = librpa_int::decide_auto_routing(n_atoms, n_kpoints * opts.nfreq);

    if (routing == LibrpaParallelRouting::RTAU)
    {
        read_Vq_full(driver_params.input_dir, "coulomb_cut_", true);
    }
    else
    {
        // NOTE: local_atpair set during read_data::read_ri.
        read_Vq_row(driver_params.input_dir, "coulomb_cut_", opts.vq_threshold, local_atpair, true);
    }
    profiler.stop("read_vq_cut");

    h.build_exx(opts);

    // TODO hide to another API
    auto pds = librpa_int::api::get_dataset_instance(h.get_c_handler());
    const auto meanfield = pds->mf;
    auto &pexx = pds->p_exx;
    pexx->build_KS_kgrid_blacs(pds->blacs_ctxt_h);

    // NOTE: since the internal object does not support parallel calculation with different k-points on each process
    // we compute all and only let master process print
    const int i_state_low = 0;
    const int n_kpoints = meanfield.get_n_kpoints();
    const int i_state_high = meanfield.get_n_bands();
    const int n_states_calc = i_state_high - i_state_low;

    std::vector<double> exx_ks(meanfield.get_n_states());
    const int n_kpoints_task = n_kpoints;
    for (int isp = 0; isp != meanfield.get_n_spins(); isp++)
    {
        for (int i_ik = 0; i_ik < n_kpoints_task; i_ik++)
        {
            // const auto ik = i_kpoints_compute[i_ik];
            const auto ik = i_ik;
            for (int ib = i_state_low; ib < i_state_high; ib++)
            {
                const int index = isp * n_kpoints_task * n_states_calc + ik * n_states_calc + ib - i_state_low;
                exx_ks[index] = pexx->Eexx[isp][ik][ib];
            }
        }
    }

    // auto exx_ks = librpa_int::app::compute_exx_orbital_energy_(i_state_low, i_state_high, -1, nullptr);

    if (exx_ks.size() > 0 && mpi_comm_global_h.is_root())
    {
        for (int isp = 0; isp != meanfield.get_n_spins(); isp++)
        {
            lib_printf("Spin channel %1d\n", isp+1);
            for (int ik = 0; ik != meanfield.get_n_kpoints(); ik++)
            {
                cout << "k-point " << ik + 1 << ": " << pds->pbc.kfrac_list[ik] << endl;
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
