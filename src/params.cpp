#include "params.h"
#include <string>
#include <utility>
#include <vector>

std::string Params::task = "rpa";
std::string Params::output_file="stdout";
std::string Params::output_dir = "librpa.d";
std::string Params::tfgrids_type = "minimax";
std::string Params::DFT_software =  "auto";
std::string Params::parallel_routing = "auto";

int Params::nfreq = 0;

double Params::gf_R_threshold = 1e-4;
double Params::cs_threshold = 1e-4;
double Params::vq_threshold = 0;
double Params::sqrt_coulomb_threshold = 1e-8;
double Params::libri_chi0_threshold_CSM = 0.0;
double Params::libri_chi0_threshold_C = 0.0;
double Params::libri_chi0_threshold_G = 0.0;
double Params::libri_exx_threshold_CSM = 0.0;
double Params::libri_exx_threshold_C = 0.0;
double Params::libri_exx_threshold_D = 0.0;
double Params::libri_exx_threshold_V = 0.0;

bool Params::binary_input = false;
bool Params::use_scalapack_ecrpa = true;
bool Params::use_scalapack_gw_wc = false;
bool Params::debug = false;
bool Params::output_gw_sigc_mat = true;
bool Params::replace_w_head = true;

int Params::option_dielect_func = 2;

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
            {"option_dielect_func", option_dielect_func},
        };

    const std::vector<std::pair<std::string, std::string>> str_params
        {
            {"task", task},
            {"output_dir", output_dir},
            {"tfgrids_type", tfgrids_type},
            {"parallel_routing", parallel_routing},
        };

    const std::vector<std::pair<std::string, bool>> bool_params
        {
            {"debug", debug},
            {"binary_input", binary_input},
            {"use_scalapack_ecrpa", use_scalapack_ecrpa},
            {"use_scalapack_gw_wc", use_scalapack_gw_wc},
            {"output_gw_sigc_mat", output_gw_sigc_mat},
            {"replace_w_head", replace_w_head},
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
