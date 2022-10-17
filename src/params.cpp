#include "params.h"
#include <string>
#include <utility>
#include <vector>

void Params::check_consistency()
{
}

void Params::print()
{
    const std::vector<std::pair<std::string, double>> double_params
        {
            {"gf_R_threshold", gf_R_threshold},
            {"cs_R_threshold", cs_threshold},
            {"vq_threshold", vq_threshold},
            {"sqrt_coulomb_threshold", sqrt_coulomb_threshold},
            {"libri_chi0_threshold_C", libri_chi0_threshold_C},
            {"libri_chi0_threshold_G", libri_chi0_threshold_G},
            {"libri_chi0_csm_threshold", libri_chi0_csm_threshold},
        };

    for (auto &param: double_params)
        printf("%s = %f\n", param.first.c_str(), param.second);
}

Params params;
