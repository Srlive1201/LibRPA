// driver headers
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

    // ofs_myid << "n_states_calc " << n_states_calc << endl;
    // ofs_myid << "n_kpoints " << n_kpoints << endl;

    pds->comm_h.barrier();
    for (int isp = 0; isp != meanfield.get_n_spins(); isp++)
    {
        std::vector<double> exx_sp_collected(n_states_calc * n_kpoints, 0.0);
        const auto exx_isp = pexx->Eexx.at(isp);
        for (const auto &[ik, exx]: exx_isp)
        {
            const int index = ik * n_states_calc;
            // ofs_myid << "ik " << ik << " index " << index << endl;
            for (int i = 0; i < n_states_calc; i++)
                exx_sp_collected[index + i] = exx.at(i_state_low + i);
        }
        if (pds->comm_h.myid == 0)
        {
            MPI_Reduce(MPI_IN_PLACE, exx_sp_collected.data(), n_states_calc * n_kpoints, MPI_DOUBLE, MPI_SUM, 0, pds->comm_h.comm);
            cout << "Spin channel " << isp+1 << endl;
            for (int ik = 0; ik < n_kpoints; ik++)
            {
                const int index = ik * n_states_calc;
                cout << "k-point " << ik + 1 << ": " << pds->pbc.kfrac_list[ik] << endl;
                lib_printf("%-4s  %-10s  %-10s\n", "Band", "e_exx (Ha)", "e_exx (eV)");
                for (int ib = 0; ib != n_states_calc; ib++)
                {
                    const auto &e = exx_sp_collected[index + ib];
                    lib_printf("%4d  %10.5f  %10.5f\n", i_state_low + ib + 1, e, HA2EV * e);
                }
                cout << endl;
            }
        }
        else
        {
            MPI_Reduce(exx_sp_collected.data(), exx_sp_collected.data(), n_states_calc * n_kpoints, MPI_DOUBLE, MPI_SUM, 0, pds->comm_h.comm);
        }
        pds->comm_h.barrier();
    }

    // auto exx_ks = librpa_int::app::compute_exx_orbital_energy_(i_state_low, i_state_high, -1, nullptr);
    profiler.stop("exx");
}
