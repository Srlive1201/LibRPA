#include "driver_utils.h"

#include <cassert>
#include <stdexcept>

#include "dielecmodel.h"
#include "driver_params.h"
#include "envs_mpi.h"
#include "fitting.h"
#include "interpolate.h"
#include "meanfield.h"
#include "read_data.h"
#include "ri.h"

std::vector<double> interpolate_dielec_func(int option, const std::vector<double> &frequencies_in,
                                            const std::vector<double> &df_in,
                                            const std::vector<double> &frequencies_target)
{
    std::vector<double> df_target;

    switch (option)
    {
        case 0: /* No extrapolation, copy the input data to target */
        {
            assert(frequencies_in.size() == frequencies_target.size());
            df_target = df_in;
            break;
        }
        case 1: /* Use spline interpolation */
        {
            df_target =
                LIBRPA::utils::interp_cubic_spline(frequencies_in, df_in, frequencies_target);
            break;
        }
        case 2: /* Use dielectric model for fitting */
        {
            LIBRPA::utils::LevMarqFitting levmarq;
            // use double-dispersion Havriliak-Negami model
            // initialize the parameters as 1.0
            std::vector<double> pars(DoubleHavriliakNegami::d_npar, 1);
            pars[0] = pars[4] = df_in[0];
            df_target =
                levmarq.fit_eval(pars, frequencies_in, df_in, DoubleHavriliakNegami::func_imfreq,
                                 DoubleHavriliakNegami::grad_imfreq, frequencies_target);
            break;
        }
        case 3: /* Read velocity matrix and calculate head and wing */
        {
            int n_basis, n_states, n_spin;
            // using LIBRPA::envs::mpi_comm_global_h;
            /*read_scf_occ_eigenvalues("./pyatb_librpa_df/band_out", pyatb_meanfield);
            read_eigenvector("./pyatb_librpa_df/", pyatb_meanfield);
            read_velocity("./pyatb_librpa_df/velocity_matrix", pyatb_meanfield);
            std::vector<Vector3_Order<double>> kfrac_band;
            kfrac_band =
                read_band_kpath_info(n_basis, n_states, n_spin, "./pyatb_librpa_df/k_path_info");

            df_headwing.set(pyatb_meanfield, kfrac_band, frequencies_target, n_basis, n_states,
                            n_spin);*/

            std::string file_abacus = driver_params.input_dir + "pyatb_librpa_df/velocity_matrix";
            std::string file_aims = driver_params.input_dir + "moment_KS_spin_01_kpt_000001.dat";
            ifstream infile_abacus;
            ifstream infile_aims;
            infile_abacus.open(file_abacus);
            infile_aims.open(file_aims);
            if (infile_abacus.is_open())
            {
                /*read_velocity(file_abacus, meanfield);
                std::vector<Vector3_Order<double>> kfrac_band = read_band_kpath_info(
                    n_basis, n_states, n_spin, "./pyatb_librpa_df/k_path_info");
                df_headwing.set(meanfield, kfrac_band, frequencies_target, n_basis, n_states,
                                n_spin);*/
                // read_scf_occ_eigenvalues("./pyatb_librpa_df/band_out", pyatb_meanfield);
                // read_eigenvector("./pyatb_librpa_df/", pyatb_meanfield);

                read_velocity(file_abacus, meanfield);
                std::vector<Vector3_Order<double>> kfrac_band =
                    read_band_kpath_info(driver_params.input_dir + "pyatb_librpa_df/k_path_info",
                                         n_basis, n_states, n_spin);
                df_headwing.set(meanfield, kfrac_band, frequencies_target, n_basis, n_states,
                                n_spin);
            }
            else if (infile_aims.is_open())
            {
                read_velocity_aims(meanfield, driver_params.input_dir);
                n_basis = meanfield.get_n_aos();
                n_states = meanfield.get_n_bands();
                n_spin = meanfield.get_n_spins();
                df_headwing.set(meanfield, kfrac_list, frequencies_target, n_basis, n_states,
                                n_spin);
            }
            else
            {
                throw std::runtime_error("Cannot find moment files for head/wing!");
            }
            infile_abacus.close();
            infile_aims.close();

            df_headwing.cal_head();
            df_headwing.cal_wing();

            df_target = df_headwing.get_head_vec();
            break;
        }
        case 4: /* Read velocity matrix and calculate head only */
        {
            int n_basis, n_states, n_spin;

            std::string file_abacus = driver_params.input_dir + "pyatb_librpa_df/velocity_matrix";
            std::string file_aims = driver_params.input_dir + "moment_KS_spin_01_kpt_000001.dat";
            ifstream infile_abacus;
            ifstream infile_aims;
            infile_abacus.open(file_abacus);
            infile_aims.open(file_aims);
            if (infile_abacus.is_open())
            {
                read_velocity(file_abacus, meanfield);
                std::vector<Vector3_Order<double>> kfrac_band = kfrac_band =
                    read_band_kpath_info(driver_params.input_dir + "pyatb_librpa_df/k_path_info",
                                         n_basis, n_states, n_spin);
                df_headwing.set(meanfield, kfrac_band, frequencies_target, n_basis, n_states,
                                n_spin);
            }
            else if (infile_aims.is_open())
            {
                read_velocity_aims(meanfield, driver_params.input_dir);
                n_basis = meanfield.get_n_aos();
                n_states = meanfield.get_n_bands();
                n_spin = meanfield.get_n_spins();
                df_headwing.set(meanfield, kfrac_list, frequencies_target, n_basis, n_states,
                                n_spin);
            }
            else
            {
                throw std::runtime_error("Cannot find moment files for head!");
            }
            infile_abacus.close();
            infile_aims.close();

            df_headwing.cal_head();

            df_target = df_headwing.get_head_vec();
            df_headwing.test_head();
            break;
        }
        default:
            throw std::logic_error("Unsupported value for option");
    }

    return df_target;
}
