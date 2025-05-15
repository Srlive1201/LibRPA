#include "app_exx.h"

#include <numeric>
#include <stdexcept>
#include <string>

#include "coulmat.h"
#include "envs_mpi.h"
#include "exx.h"
#include "meanfield.h"
#include "params.h"
#include "pbc.h"
#include "ri.h"

namespace LIBRPA
{

namespace app
{

std::vector<double> compute_exx_orbital_energy_(int i_state_low, int i_state_high,
                                                int n_kpoints_task, const int *i_kpoints_task)
{
    using LIBRPA::envs::mpi_comm_global_h;
    using LIBRPA::utils::lib_printf;

    std::vector<int> i_kpoints_compute(0);

    // Determine the indices of k-points to compute
    if (n_kpoints_task < 0)
    {
        // Compute all k-points
        n_kpoints_task = meanfield.get_n_kpoints();
        i_kpoints_compute.resize(n_kpoints_task);
        std::iota(i_kpoints_compute.begin(), i_kpoints_compute.end(), 0);
    }
    else
    {
        // Compute only specified k-points
        i_kpoints_compute = std::vector<int>(i_kpoints_task, i_kpoints_task + n_kpoints_task);
    }

    const auto n_spins = meanfield.get_n_spins();

    // Initialize the return vector with correct size
    std::vector<double> exx_state(0);
    const int n_states_calc = i_state_high - i_state_low;
    if (n_states_calc > 0)
    {
        exx_state.resize(n_spins * n_kpoints_task * (i_state_high - i_state_low));
    }
    else if (i_state_low >= meanfield.get_n_bands())
    {
        throw std::runtime_error(
            std::string(__FILE__) + ":" + std::to_string(__LINE__) + ":" +
            std::string(__FUNCTION__) + ": " + "i_state_low should be no larger than n_bands = " +
            std::to_string(meanfield.get_n_bands()) + ", got " + std::to_string(i_state_low));
    }

    Vector3_Order<int> period{kv_nmp[0], kv_nmp[1], kv_nmp[2]};
    auto Rlist = construct_R_grid(period);

    const auto VR = FT_Vq(Vq_cut, meanfield.get_n_kpoints(), Rlist, true);
    // TODO: kfrac_list should depend on i_kpoints_compute
    auto exx = LIBRPA::Exx(meanfield, kfrac_list, period);
    if (Params::use_shrink_abfs)
    {
        if (Params::use_soc)
            exx.build<std::complex<double>>(Cs_shrinked_data, Rlist, VR);
        else
            exx.build<double>(Cs_shrinked_data, Rlist, VR);
    }
    else
    {
        if (Params::use_soc)
            exx.build<std::complex<double>>(Cs_data, Rlist, VR);
        else
            exx.build<double>(Cs_data, Rlist, VR);
    }
    exx.build_KS_kgrid();

    for (int isp = 0; isp != meanfield.get_n_spins(); isp++)
    {
        for (int i_ik = 0; i_ik < n_kpoints_task; i_ik++)
        {
            const auto ik = i_kpoints_compute[i_ik];
            for (int ib = i_state_low; ib < i_state_high; ib++)
            {
                const int index =
                    isp * n_kpoints_task * n_states_calc + ik * n_states_calc + ib - i_state_low;
                exx_state[index] = exx.Eexx[isp][ik][ib];
            }
        }
    }
    return exx_state;
}

} /* end of namespace app */

} /* end of namespace LIBRPA */
