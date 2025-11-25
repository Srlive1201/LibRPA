// Public headers (prefixed by librpa)
#include "librpa_options.h"
#include "librpa_enums.h"

#include <string>
#include <cstring>
#include <stdexcept>

// C APIs
void librpa_init_options(LibrpaOptions *opts)
{
    librpa_set_output_dir(opts, ".");

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

    opts->output_gw_sigc_mat = LIBRPA_SWITCH_OFF;
    opts->output_gw_sigc_mat_rt = LIBRPA_SWITCH_OFF;
    opts->output_gw_sigc_mat_rf = LIBRPA_SWITCH_OFF;
}

static std::string path_as_directory(const std::string &dirpath)
{
    if (dirpath.find(":") != std::string::npos)
    {
        throw std::runtime_error("dirpath contains invalid character (:) for POSIX path");
    }

    if (dirpath.back() != '/')
    {
        return dirpath + '/';
    }

    return dirpath;
}

void librpa_set_output_dir(LibrpaOptions *opts, const char *output_dir)
{
    std::string output_dir_s = path_as_directory(output_dir);
    strcpy(opts->output_dir, output_dir_s.c_str());
}
