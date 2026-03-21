#include "driver.h"

#include <ios>
#include <sstream>
#include <string>
#include "librpa_enums.h"

namespace driver
{

DriverParams::DriverParams():
    task("rpa"),
    constants_choice("internal"),
    input_dir(""),
    output_gw_spec_func(false),
    sf_omega_start(0.0),
    sf_omega_end(1.0),
    sf_omega_step(0.1),
    sf_gf_omega_shift(0.01),
    sf_sigc_omega_shift(0.01),
    sf_state_start(0),
    sf_state_end(10000)
{
}

std::string DriverParams::format()
{
    std::stringstream ss;
    ss << "task = " << task << std::endl;
    ss << "input_dir = " << input_dir << std::endl;
    ss << "constants_choice = " << constants_choice << std::endl;
    ss << "output_gw_spec_func = " << std::boolalpha << output_gw_spec_func << std::endl;
    if (output_gw_spec_func)
    {
        ss << "sf_omega_start = " << sf_omega_start << std::endl;
        ss << "sf_omega_end   = " << sf_omega_end << std::endl;
        ss << "sf_omega_step  = " << sf_omega_step << std::endl;
    }
    return ss.str();
}

DriverParams driver_params;

const std::string input_filename = "librpa.in";

std::vector<int> atom_types;
size_t n_atoms;

int n_spins = 0;
int n_kpoints = 0;
int n_kpoints_band = 0;
int n_ibz_kpoints = 0;
int n_states = 0;
int n_basis_wfc = 0;

std::vector<int> iks_eigvec_this;
std::vector<int> iks_band_eigvec_this;

// Used to read Coulomb matrix data.
// Should be consistent with the internal `atpairs_local` of the Dataset object
std::vector<std::pair<size_t, size_t>> local_atpair;

std::vector<librpa_int::Vector3_Order<double>> ibz_kpoints;
std::vector<librpa_int::Vector3_Order<double>> kfrac_band;

librpa::Handler h;

librpa::Options opts;

#define normal_pair(name) {#name, opts.name}
#define bool_pair(name) {#name, get_bool(opts.name)}

std::string format_runtime_options(const librpa::Options &opts) noexcept
{
    std::stringstream ss;
    const std::vector<std::pair<std::string, double>> double_params
        {
            {"gf_R_threshold", opts.gf_threshold},
            {"cs_R_threshold", opts.cs_threshold},
            normal_pair(vq_threshold),
            normal_pair(sqrt_coulomb_threshold),
            normal_pair(libri_chi0_threshold_C),
            normal_pair(libri_chi0_threshold_G),
            normal_pair(libri_exx_threshold_C),
            normal_pair(libri_exx_threshold_D),
            normal_pair(libri_exx_threshold_V),
            normal_pair(libri_g0w0_threshold_C),
            normal_pair(libri_g0w0_threshold_G),
            normal_pair(libri_g0w0_threshold_Wc),
        };

    const std::vector<std::pair<std::string, int>> int_params
        {
            normal_pair(nfreq),
            normal_pair(n_params_anacon),
            normal_pair(option_dielect_func),
        };

    const std::vector<std::pair<std::string, std::string>> str_params
        {
            {"output_dir", opts.output_dir},
            {"output_level", get_verbose_string(opts.output_level)},
            {"tfgrids_type", get_tfgrid_string(opts.tfgrids_type)},
            {"parallel_routing", get_routing_string(opts.parallel_routing)},
        };

    const std::vector<std::pair<std::string, bool>> bool_params
        {
            {"debug", opts.output_level >= LIBRPA_VERBOSE_DEBUG},
            bool_pair(replace_w_head),
            bool_pair(use_scalapack_ecrpa),
            bool_pair(use_scalapack_gw_wc),
            bool_pair(use_kpara_scf_eigvec),
            bool_pair(output_gw_sigc_mat),
            bool_pair(output_gw_sigc_mat_rf),
            bool_pair(output_gw_sigc_mat_rt),
        };

    for (const auto &param: str_params)
        ss << param.first.c_str() << " = " << param.second << std::endl;

    for (const auto &param: int_params)
        ss << param.first.c_str() << " = " << param.second << std::endl;

    for (const auto &param: double_params)
        ss << param.first.c_str() << " = " << param.second << std::endl;

    for (const auto &param: bool_params)
        ss << param.first.c_str() << " = " << std::boolalpha << param.second << std::endl;

    return ss.str();
}

#undef normal_pair
#undef bool_pair

LibrpaTimeFreqGrid get_tfgrid_type(const std::string& grid_str)
{
    if (grid_str == "unset")
        return LibrpaTimeFreqGrid::TFGRID_UNSET;
    if (grid_str == "GL")
        return LibrpaTimeFreqGrid::GaussLegendre;
    if (grid_str == "GC-I")
        return LibrpaTimeFreqGrid::GaussChebyshevI;
    if (grid_str == "GL-II")
        return LibrpaTimeFreqGrid::GaussChebyshevII;
    if (grid_str == "minimax")
        return LibrpaTimeFreqGrid::Minimax;
    if (grid_str == "evenspaced")
        return LibrpaTimeFreqGrid::EvenSpaced;
    if (grid_str == "evenspaced_tf")
        return LibrpaTimeFreqGrid::EvenSpaced_TF;
    throw std::runtime_error("Unknown time-frequency grid string: " + grid_str);
}

std::string get_tfgrid_string(const LibrpaTimeFreqGrid& grid_type) noexcept
{
    if (grid_type == LibrpaTimeFreqGrid::GaussLegendre)
        return "GL";
    if (grid_type == LibrpaTimeFreqGrid::GaussChebyshevI)
        return "GC-I";
    if (grid_type == LibrpaTimeFreqGrid::GaussChebyshevII)
        return "GL-II";
    if (grid_type == LibrpaTimeFreqGrid::Minimax)
        return "minimax";
    if (grid_type == LibrpaTimeFreqGrid::EvenSpaced)
        return "evenspaced";
    if (grid_type == LibrpaTimeFreqGrid::EvenSpaced_TF)
        return "evenspaced_tf";
    return "unset";
}

LibrpaParallelRouting get_parallel_routing(const std::string& routing_str_low)
{
    if (routing_str_low == "auto") return LibrpaParallelRouting::AUTO;
    if (routing_str_low == "rtau") return LibrpaParallelRouting::RTAU;
    if (routing_str_low == "atompair") return LibrpaParallelRouting::ATOMPAIR;
    if (routing_str_low == "libri") return LibrpaParallelRouting::LIBRI;
    throw std::runtime_error("Unknown parallel routing string: " + routing_str_low);
}

std::string get_routing_string(LibrpaParallelRouting routing)
{
    if (routing == LibrpaParallelRouting::AUTO) return "auto";
    if (routing == LibrpaParallelRouting::ATOMPAIR) return "atompair";
    if (routing == LibrpaParallelRouting::RTAU) return "rtau";
    if (routing == LibrpaParallelRouting::LIBRI) return "libri";
    return "unset";
}

std::string get_verbose_string(LibrpaVerbose verbose)
{
    if (verbose == LIBRPA_VERBOSE_DEBUG) return "debug";
    if (verbose == LIBRPA_VERBOSE_WARN) return "warn";
    if (verbose == LIBRPA_VERBOSE_INFO) return "info";
    if (verbose == LIBRPA_VERBOSE_CRITICAL) return "critical";
    return "silent";
}

LibrpaVerbose get_verbose(const std::string& verbose_str_low)
{
    if (verbose_str_low == "debug") return LIBRPA_VERBOSE_DEBUG;
    if (verbose_str_low == "warn") return LIBRPA_VERBOSE_WARN;
    if (verbose_str_low == "warning") return LIBRPA_VERBOSE_WARN;
    if (verbose_str_low == "info") return LIBRPA_VERBOSE_INFO;
    if (verbose_str_low == "critical") return LIBRPA_VERBOSE_CRITICAL;
    if (verbose_str_low == "silent") return LIBRPA_VERBOSE_SILENT;
    throw std::runtime_error("Unknown verbose level string: " + verbose_str_low);
}

}
