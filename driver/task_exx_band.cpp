#include "task_exx_band.h"

#include <fstream>
#include <sstream>

#include "../src/core/coulmat.h"
#include "../src/core/exx.h"
#include "../src/core/meanfield.h"
#include "../src/core/pbc.h"
#include "../src/core/ri.h"
#include "../src/mpi/envs_mpi.h"
#include "../src/utils/constants.h"
#include "../src/utils/envs_io.h"
#include "../src/utils/profiler.h"
#include "driver_params.h"
#include "read_data.h"

void task_exx_band()
{
    using LIBRPA::envs::mpi_comm_global_h;
    using LIBRPA::envs::ofs_myid;
    using LIBRPA::utils::lib_printf;

    Profiler::start("exx_band");

    Vector3_Order<int> period {kv_nmp[0], kv_nmp[1], kv_nmp[2]};
    auto Rlist = construct_R_grid(period);

    vector<Vector3_Order<double>> qlist;
    for ( auto q_weight: irk_weight)
    {
        qlist.push_back(q_weight.first);
    }

    // Prepare time-frequency grids
    Profiler::start("read_vq_cut", "Load truncated Coulomb");
    read_Vq_full(driver_params.input_dir, "coulomb_cut_", true);
    Profiler::stop("read_vq_cut");

    std::vector<double> epsmac_LF_imagfreq_re;

    Profiler::start("exx_real_space", "Build exchange self-energy");
    auto exx = LIBRPA::Exx(meanfield, kfrac_list, period);
    {
        Profiler::start("ft_vq_cut", "Fourier transform truncated Coulomb");
        const auto VR = FT_Vq(Vq_cut, meanfield.get_n_kpoints(), Rlist, true);
        Profiler::stop("ft_vq_cut");

        Profiler::start("exx_real_work");
        exx.build(Cs_data, Rlist, VR);
        Profiler::stop("exx_real_work");
    }
    Profiler::stop("exx_real_space");
    std::flush(ofs_myid);

    Profiler::start("read_vxc", "Load DFT xc potential");
    std::vector<matrix> vxc;
    int flag_read_vxc = read_vxc(driver_params.input_dir + "vxc_out", vxc);
    Profiler::stop("read_vxc");

    if (flag_read_vxc != 0)
    {
        if (mpi_comm_global_h.myid == 0)
        {
            lib_printf("Error in reading Vxc on kgrid, task failed!\n");
        }
        Profiler::stop("exx_band");
        return;
    }

    Profiler::start("exx_ks_kgrid");
    exx.build_KS_kgrid();
    Profiler::stop("exx_ks_kgrid");

    if (mpi_comm_global_h.is_root())
    {
        // display results
        const std::string banner(90, '-');
        printf("Printing EXX@DFT energy [unit: eV]\n\n");
        for (int i_spin = 0; i_spin < meanfield.get_n_spins(); i_spin++)
        {
            for (int i_kpoint = 0; i_kpoint < meanfield.get_n_kpoints(); i_kpoint++)
            {
                const auto &k = kfrac_list[i_kpoint];
                printf("spin %2d, k-point %4d: (%.5f, %.5f, %.5f) \n",
                        i_spin+1, i_kpoint+1, k.x, k.y, k.z);
                printf("%90s\n", banner.c_str());
                printf("%5s %16s %16s %16s %16s %16s\n", "State", "occ", "e_mf", "v_xc", "v_exx", "e_exx");
                printf("%90s\n", banner.c_str());
                for (int i_state = 0; i_state < meanfield.get_n_bands(); i_state++)
                {
                    const auto &occ_state = meanfield.get_weight()[i_spin](i_kpoint, i_state) * meanfield.get_n_kpoints();
                    const auto &eks_state = meanfield.get_eigenvals()[i_spin](i_kpoint, i_state) * HA2EV;
                    const auto &exx_state = exx.Eexx[i_spin][i_kpoint][i_state] * HA2EV;
                    const auto &vxc_state = vxc[i_spin](i_kpoint, i_state) * HA2EV;
                    printf("%5d %16.5f %16.5f %16.5f %16.5f %16.5f\n",
                           i_state+1, occ_state, eks_state, vxc_state, exx_state, eks_state - vxc_state + exx_state);
                }
                printf("\n");
            }
        }
    }

    /*
     * Compute the EXX energies on band k-paths
     */
    // Reset k-space EXX and Sigmac matrices to avoid warning from internal reset
    exx.reset_kspace();

    /* Below we handle the band k-points data
     * First load the information of k-points along the k-path */
    Profiler::start("exx_band_load_band_mf");
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

    /* Set the same Fermi energy as in SCF */
    meanfield_band.get_efermi() = meanfield.get_efermi();
    Profiler::stop("exx_band_load_band_mf");

    Profiler::start("exx_rotate_KS");
    exx.build_KS_band(meanfield_band.get_eigenvectors(), kfrac_band);
    Profiler::stop("exx_rotate_KS");
    std::flush(ofs_myid);

    Profiler::start("read_vxc", "Load DFT xc potential");
    auto vxc_band = read_vxc_band(driver_params.input_dir, n_states_band, n_spin_band, kfrac_band.size());
    Profiler::stop("read_vxc");
    std::flush(ofs_myid);

    // TODO: parallelize analytic continuation and QPE solver among tasks
    if (mpi_comm_global_h.is_root())
    {
        const auto &mf = meanfield_band;
        // display results
        for (int i_spin = 0; i_spin < mf.get_n_spins(); i_spin++)
        {
            std::ofstream ofs_ks;
            std::ofstream ofs_hf;
            std::stringstream fn;

            fn << "EXX_band_spin_" << i_spin + 1 << ".dat";
            ofs_hf.open(fn.str());

            fn.str("");
            fn.clear();
            fn << "KS_band_spin_" << i_spin + 1 << ".dat";
            ofs_ks.open(fn.str());

            ofs_hf << std::fixed;
            ofs_ks << std::fixed;

            for (int i_kpoint = 0; i_kpoint < mf.get_n_kpoints(); i_kpoint++)
            {
                const auto &k = kfrac_band[i_kpoint];
                ofs_ks << std::setw(5) << i_kpoint + 1 << std::setw(15) << std::setprecision(7) << k.x << std::setw(15) << std::setprecision(7) << k.y << std::setw(15) << std::setprecision(7) << k.z;
                ofs_hf << std::setw(5) << i_kpoint + 1 << std::setw(15) << std::setprecision(7) << k.x << std::setw(15) << std::setprecision(7) << k.y << std::setw(15) << std::setprecision(7) << k.z;
                for (int i_state = 0; i_state < meanfield.get_n_bands(); i_state++)
                {
                    const auto &occ_state = mf.get_weight()[i_spin](i_kpoint, i_state);
                    const auto &eks_state = mf.get_eigenvals()[i_spin](i_kpoint, i_state) * HA2EV;
                    const auto &exx_state = exx.Eexx[i_spin][i_kpoint][i_state] * HA2EV;
                    const auto &vxc_state = vxc_band[i_spin](i_kpoint, i_state) * HA2EV;
                    // const auto &resigc = sigc_all[i_spin][i_kpoint][i_state].real() * HA2EV;
                    // const auto &imsigc = sigc_all[i_spin][i_kpoint][i_state].imag() * HA2EV;
                    ofs_ks << std::setw(15) << std::setprecision(5) << occ_state << std::setw(15) << std::setprecision(5) << eks_state;
                    ofs_hf << std::setw(15) << std::setprecision(5) << occ_state << std::setw(15) << std::setprecision(5) << eks_state - vxc_state + exx_state;
                }
                ofs_hf << "\n";
                ofs_ks << "\n";
            }
        }
    }


    Profiler::stop("exx_band");
}
