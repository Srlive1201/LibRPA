// Internal headers
#include "dataset_helper.h"
#include "dataset.h"

namespace librpa_int
{

void initialize_tfgrids(Dataset &ds, const LibrpaOptions &opts)
{
    ds.tfg.reset(opts.nfreq);
    double emin = opts.tfgrids_freq_min;
    double eintv = opts.tfgrids_freq_interval;
    double emax = opts.tfgrids_freq_max;
    double tmin = opts.tfgrids_time_min;
    double tintv = opts.tfgrids_time_interval;
    if (opts.tfgrids_type == LibrpaTimeFreqGrid::Minimax)
    {
        ds.mf.get_E_min_max(emin, emax);
    }
    ds.tfg.generate(opts.tfgrids_type, emin, eintv, emax, tmin, tintv);
}

}
