// Public headers (prefixed by librpa)
#include "librpa_enums.h"
#include "librpa_options.h"

#include "../io/fs.h"

#include <string>
#include <cstring>

// C APIs
void librpa_init_options(LibrpaOptions *opts)
{
    librpa_set_output_dir(opts, ".");

    opts->parallel_routing = LibrpaParallelRouting::AUTO;
    opts->output_level = LIBRPA_VERBOSE_INFO;
    opts->vq_threshold = 0.0e0;
    opts->use_soc = LIBRPA_SWITCH_OFF;
    opts->use_kpara_scf_eigvec = LIBRPA_SWITCH_OFF;

    opts->tfgrids_type = LibrpaTimeFreqGrid::TFGRID_UNSET;
    opts->nfreq = 6;
    opts->tfgrids_freq_min = 0.005;
    opts->tfgrids_freq_interval = 0.0;
    opts->tfgrids_freq_max = 1000.0;
    opts->tfgrids_time_min = 0.005;
    opts->tfgrids_time_interval = 0.0;

    opts->gf_threshold = 0.0e0;
    opts->use_scalapack_ecrpa = LIBRPA_SWITCH_OFF;

    opts->n_params_anacon = -1;
    opts->use_scalapack_gw_wc = LIBRPA_SWITCH_ON;
    opts->replace_w_head = LIBRPA_SWITCH_OFF;
    opts->option_dielect_func = 0;
    opts->sqrt_coulomb_threshold = 0.0e0;

    opts->libri_chi0_threshold_C = 0.0e0;
    opts->libri_chi0_threshold_G = 0.0e0;
    opts->libri_exx_threshold_C = 0.0e0;
    opts->libri_exx_threshold_D = 0.0e0;
    opts->libri_exx_threshold_V = 0.0e0;
    opts->libri_g0w0_threshold_C = 0.0e0;
    opts->libri_g0w0_threshold_G = 0.0e0;
    opts->libri_g0w0_threshold_Wc = 0.0e0;

    opts->output_gw_sigc_mat = LIBRPA_SWITCH_OFF;
    opts->output_gw_sigc_mat_rt = LIBRPA_SWITCH_OFF;
    opts->output_gw_sigc_mat_rf = LIBRPA_SWITCH_OFF;
}

void librpa_set_output_dir(LibrpaOptions *opts, const char *output_dir)
{
    std::string output_dir_s = librpa_int::path_as_directory(output_dir);
    strcpy(opts->output_dir, output_dir_s.c_str());
}
