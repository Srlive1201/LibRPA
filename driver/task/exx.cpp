// driver headers
#include <cstring>
#include <iostream>
#include <ostream>
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
    using std::cout;
    using std::endl;
    using std::map;
    using namespace librpa_int::global;
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
    const int i_state_low = 0;
    const int i_state_high = n_states;
    const int n_states_calc = i_state_high - i_state_low;
    const auto exx_ks = h.get_exx_pot_kgrid(opts, n_spins, iks_eigvec_this, i_state_low, i_state_high);

    const auto &kfrac_list = librpa_int::api::get_dataset_instance(h)->pbc.kfrac_list;

    mpi_comm_global_h.barrier();
    for (int isp = 0; isp != n_spins; isp++)
    {
        std::vector<double> exx_sp_collected(n_states_calc * n_kpoints, 0.0);
        const int st = isp * n_states_calc * iks_eigvec_this.size();
        for (size_t ik_this = 0; ik_this < iks_eigvec_this.size(); ik_this++)
        {
            const int ik = iks_eigvec_this[ik_this];
            const int index_collect = ik * n_states_calc;
            const int index = st + ik_this * n_states_calc;
            memcpy(exx_sp_collected.data() + index_collect, exx_ks.data() + index, n_states_calc * sizeof(double));
        }
        mpi_comm_global_h.reduce(MPI_IN_PLACE, exx_sp_collected.data(), n_states_calc * n_kpoints, 0, MPI_SUM);
        if (mpi_comm_global_h.myid == 0)
        {
            cout << "Spin channel " << isp+1 << endl;
            for (int ik = 0; ik < n_kpoints; ik++)
            {
                const int index = ik * n_states_calc;
                cout << "k-point " << ik + 1 << ": " << kfrac_list[ik] << endl;
                lib_printf("%-4s  %-10s  %-10s\n", "Band", "e_exx (Ha)", "e_exx (eV)");
                for (int ib = 0; ib != n_states_calc; ib++)
                {
                    const auto &e = exx_sp_collected[index + ib];
                    lib_printf("%4d  %10.5f  %10.5f\n", i_state_low + ib + 1, e, HA2EV * e);
                }
                cout << endl;
            }
        }
        mpi_comm_global_h.barrier();
    }

    // auto exx_ks = librpa_int::app::compute_exx_orbital_energy_(i_state_low, i_state_high, -1, nullptr);
    profiler.stop("exx");
}
