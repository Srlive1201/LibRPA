#include "utils_timefreq.h"


namespace LIBRPA
{

namespace utils
{

TFGrids generate_timefreq_grids(unsigned ngrids, const std::string &grid_type_str, const MeanField &mf)
{
    TFGrids tfg(ngrids);
    auto grid_type = TFGrids::get_grid_type(grid_type_str);

    double emin = -1;
    double eintv = -1;
    double tmin = -1;
    double tintv = -1;
    double emax = -1;

    switch (grid_type)
    {
        case (TFGrids::GRID_TYPES::Minimax):
        {
            mf.get_E_min_max(emin, emax);
            break;
        }
        case (TFGrids::GRID_TYPES::EvenSpaced):
        {
            emin = 0.005;
            eintv = 0.0;
            break;
        }
        case (TFGrids::GRID_TYPES::EvenSpaced_TF):
        {
            emin = 0.005;
            eintv = 0.0;
            tmin = 0.005;
            tintv = 0.0;
            break;
        }
        default:
        // Other cases are handled within TFGrids
        {}
    }

    tfg.generate(grid_type, emin, eintv, emax, tmin, tintv);

    return tfg;
}

} /* end of namespace utils */

} /* end of namespace LIBRPA */
