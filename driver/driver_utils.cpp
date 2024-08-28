#include "driver_utils.h"

#include <cassert>
#include <stdexcept>

#include "fitting.h"
#include "dielecmodel.h"
#include "interpolate.h"


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
            df_target = LIBRPA::utils::interp_cubic_spline(
                    frequencies_in, df_in, frequencies_target);
            break;
        }
        case 2: /* Use dielectric model for fitting */
        {
            LIBRPA::utils::LevMarqFitting levmarq;
            // use double-dispersion Havriliak-Negami model
            // initialize the parameters as 1.0
            std::vector<double> pars(DoubleHavriliakNegami::d_npar, 1);
            pars[0] = pars[4] = df_in[0];
            df_target = levmarq.fit_eval(pars, frequencies_in, df_in,
                                         DoubleHavriliakNegami::func_imfreq,
                                         DoubleHavriliakNegami::grad_imfreq,
                                         frequencies_target);
            break;
        }
        default:
            throw std::logic_error("Unsupported value for option");
    }

    return df_target;
}
