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
            {"libri_chi0_threshold_CSM", libri_chi0_threshold_CSM},
            {"libri_exx_threshold_C", libri_exx_threshold_C},
            {"libri_exx_threshold_D", libri_exx_threshold_D},
            {"libri_exx_threshold_V", libri_exx_threshold_V},
            {"libri_exx_threshold_CSM", libri_exx_threshold_CSM},
        };

    const std::vector<std::pair<std::string, int>> int_params
        {
            {"nfreq", nfreq},
        };

    const std::vector<std::pair<std::string, std::string>> str_params
        {
            {"task", task},
            {"tfgrids_type", tfgrids_type},
            {"chi_parallel_routing", chi_parallel_routing},
            {"exx_parallel_routing", exx_parallel_routing},
        };

    const std::vector<std::pair<std::string, bool>> bool_params
        {
            {"use_libri_chi0", use_libri_chi0},
            {"use_libri_exx", use_libri_exx},
            {"use_scalapack_ecrpa", use_scalapack_ecrpa},
            {"use_scalapack_gw_wc", use_scalapack_gw_wc},
        };

    for (const auto &param: str_params)
        printf("%s = %s\n", param.first.c_str(), param.second.c_str());

    for (const auto &param: int_params)
        printf("%s = %d\n", param.first.c_str(), param.second);

    for (const auto &param: double_params)
        printf("%s = %f\n", param.first.c_str(), param.second);

    for (const auto &param: bool_params)
        printf("%s = %s\n", param.first.c_str(), param.second? "T": "F");
}

Params params;
