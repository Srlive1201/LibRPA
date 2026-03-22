#include <fstream>
#include <sstream>

#include "../task.h"
#include "../driver.h"
#include "../read_data.h"

#include "../../src/utils/profiler.h"
#include "../../src/utils/constants.h"
#include "../../src/io/global_io.h"

#include "../../src/api/instance_manager.h"

void driver::task_exx_band()
{
    using std::cout;
    using std::endl;
    using std::setw;
    using std::map;
    using librpa_int::HA2EV;
    using namespace librpa_int::global;

    profiler.start("exx_band", "Exact-exchange-only calculation of band structure");

    profiler.start("read_vq_cut", "Load truncated Coulomb");
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

    profiler.start("read_vxc", "Load DFT xc potential");
    std::vector<matrix> vxc;
    int flag_read_vxc = read_vxc(driver_params.input_dir + "vxc_out", vxc);
    profiler.stop("read_vxc");
    if (flag_read_vxc != 0)
    {
        if (mpi_comm_global_h.myid == 0)
        {
            lib_printf("Error in reading Vxc on kgrid, task failed!\n");
        }
        profiler.stop("exx_band");
        return;
    }

    h.build_exx(opts);

    const int i_state_low = 0;
    const int i_state_high = n_states;
    const int n_states_calc = i_state_high - i_state_low;

    const auto exx_ks = h.get_exx_pot_kgrid(opts, n_spins, iks_eigvec_this, i_state_low, i_state_high);

    mpi_comm_global_h.barrier();
    const std::string banner(90, '-');
    if (mpi_comm_global_h.is_root())
    {
        cout << banner << endl;
        cout << "Printing orbital energy [unit: eV]" << endl << endl;
    }
    mpi_comm_global_h.barrier();

    // Access meanfield data from the database
    const auto &mf = librpa_int::api::get_dataset_instance(h)->mf;
    const auto &kfrac_list = librpa_int::api::get_dataset_instance(h)->pbc.kfrac_list;

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
            for (int ik = 0; ik < n_kpoints; ik++)
            {
                const auto &k = kfrac_list[ik];
                cout << "Spin channel " << setw(2) << isp + 1 << ", k-point " << setw(4) << ik + 1
                     << ": (" << std::fixed << std::setprecision(5) << k.x << ", " << k.y << ", " << k.z << ")" << endl;
                cout << banner << endl;
                cout << setw(5) << "State" << setw(16) << "occ" << setw(16) << "e_mf"
                     << setw(16) << "v_xc" << setw(16) << "v_exx" << setw(16) << "e_exx" << endl;
                cout << banner << endl;
                const int index = ik * n_states_calc;
                for (int ib = 0; ib != n_states_calc; ib++)
                {
                    const int i_state = ib + i_state_low;
                    const auto &occ_state = mf.get_weight()[isp](ik, i_state) * mf.get_n_kpoints();
                    const auto &eks_state = mf.get_eigenvals()[isp](ik, i_state) * HA2EV;
                    const auto &exx_state = exx_sp_collected[index + ib] * HA2EV;
                    const auto &vxc_state = vxc[isp](ik, i_state) * HA2EV;
                    printf("%5d %16.5f %16.5f %16.5f %16.5f %16.5f\n",
                           i_state+1, occ_state, eks_state, vxc_state, exx_state, eks_state - vxc_state + exx_state);
                }
                cout << endl;
            }
        }
        mpi_comm_global_h.barrier();
    }

    /* Below we handle the band k-points data
     * First load the information of k-points along the k-path */
    profiler.start("exx_band_load_band_mf");
    read_band_kpath_info(driver_params.input_dir + "band_kpath_info");
    const int nkpts_band = kfrac_band.size();
    if (mpi_comm_global_h.is_root())
    {
        std::cout << "Band k-points to compute:\n";
        for (int ik = 0; ik < nkpts_band; ik++)
        {
            const auto &k = kfrac_band[ik];
            lib_printf("%5d %12.7f %12.7f %12.7f\n", ik + 1, k.x, k.y, k.z);
        }
    }
    mpi_comm_global_h.barrier();
    read_band_meanfield_data(driver_params.input_dir);
    profiler.stop("exx_band_load_band_mf");

    profiler.start("read_vxc_band", "Load DFT xc potential");
    const auto vxc_band = read_vxc_band(driver_params.input_dir, n_states, n_spins, kfrac_band.size());
    profiler.stop("read_vxc_band");

    const auto exx_ks_band = h.get_exx_pot_band_k(opts, n_spins, iks_band_eigvec_this, i_state_low, i_state_high);

    profiler.start("output_exx_band");
    const auto &mf_band = librpa_int::api::get_dataset_instance(h)->mf_band;
    for (int isp = 0; isp != n_spins; isp++)
    {
        std::vector<double> exx_sp_collected(n_states_calc * n_kpoints_band, 0.0);
        const int st = isp * n_states_calc * iks_band_eigvec_this.size();
        for (size_t ik_this = 0; ik_this < iks_band_eigvec_this.size(); ik_this++)
        {
            const int ik = iks_band_eigvec_this[ik_this];
            const int index_collect = ik * n_states_calc;
            const int index = st + ik_this * n_states_calc;
            memcpy(exx_sp_collected.data() + index_collect, exx_ks.data() + index, n_states_calc * sizeof(double));
        }
        mpi_comm_global_h.reduce(MPI_IN_PLACE, exx_sp_collected.data(), n_states_calc * n_kpoints_band, 0, MPI_SUM);

        if (mpi_comm_global_h.myid == 0)
        {
            std::ofstream ofs_ks, ofs_hf;
            std::stringstream fn_ks, fn_hf;
            fn_ks << "KS_band_spin_" << isp + 1 << ".dat";
            fn_hf << "EXX_band_spin_" << isp + 1 << ".dat";
            ofs_hf.open(fn_hf.str());
            ofs_ks.open(fn_ks.str());
            ofs_hf << std::fixed;
            ofs_ks << std::fixed;

            for (int ik = 0; ik < n_kpoints_band; ik++)
            {
                const int index_k = ik * n_states_calc;
                const auto &k = kfrac_band[ik];
                ofs_ks << std::setw(5) << ik + 1 << std::setw(15) << std::setprecision(7) << k.x << std::setw(15) << std::setprecision(7) << k.y << std::setw(15) << std::setprecision(7) << k.z;
                ofs_hf << std::setw(5) << ik + 1 << std::setw(15) << std::setprecision(7) << k.x << std::setw(15) << std::setprecision(7) << k.y << std::setw(15) << std::setprecision(7) << k.z;
                for (int ib = 0; ib < n_states_calc; ib++)
                {
                    const int i_state = i_state_low + ib;
                    const auto &occ_state = mf_band.get_weight()[isp](ik, i_state) * n_kpoints_band;
                    const auto &eks_state = mf_band.get_eigenvals()[isp](ik, i_state) * HA2EV;
                    const auto &exx_state = exx_sp_collected[index_k+ib] * HA2EV;
                    const auto &vxc_state = vxc_band[isp](ik, i_state) * HA2EV;
                    ofs_ks << std::setw(15) << std::setprecision(5) << occ_state << std::setw(15) << std::setprecision(5) << eks_state;
                    ofs_hf << std::setw(15) << std::setprecision(5) << occ_state << std::setw(15) << std::setprecision(5) << eks_state - vxc_state + exx_state;
                }
                ofs_hf << "\n";
                ofs_ks << "\n";
            }

            ofs_hf.close();
            ofs_ks.close();
        }
        mpi_comm_global_h.barrier();
    }
    profiler.stop("output_exx_band");

    profiler.stop("exx_band");
}
