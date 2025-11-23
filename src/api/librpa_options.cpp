// Public headers (prefixed by librpa)
#include "../../include/librpa_options.h"

// C APIs
void librpa_init_options(LibrpaOptions *opts)
{
    opts->parallel_routing = LibrpaParallelRouting::ROUTING_UNSET;
    opts->output_level = LIBRPA_VERBOSE_INFO;
    opts->tfgrids_type = LibrpaTimeFreqGrid::TFGRID_UNSET;
    opts->nfreq = 6;
    opts->cs_threshold = 0.0e0;
    opts->vq_threshold = 0.0e0;
    opts->use_soc = LIBRPA_SWITCH_OFF;

    opts->gf_threshold = 0.0e0;
    opts->use_scalapack_ecrpa = LIBRPA_SWITCH_OFF;

    opts->n_params_anacon = 6;
    opts->use_scalapack_gw_wc = LIBRPA_SWITCH_OFF;
    opts->replace_w_head = LIBRPA_SWITCH_OFF;
    opts->option_dielect_func = 0;
    opts->sqrt_coulomb_threshold = 0.0e0;

    opts->libri_chi0_threshold_C = 0.0e0;
    opts->libri_chi0_threshold_G = 0.0e0;
    opts->libri_exx_threshold_C = 0.0e0;
    opts->libri_exx_threshold_D = 0.0e0;
    opts->libri_exx_threshold_V = 0.0e0;
    opts->libri_gw_threshold_C = 0.0e0;
    opts->libri_gw_threshold_G = 0.0e0;
    opts->libri_gw_threshold_Wc = 0.0e0;
}
