#include "inputfile.h"

#include <regex>
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>

#include "driver.h"
#include "librpa_enums.h"

static const std::string SPACE_SEP = "[ \r\f\t]*";

const std::string InputParser::KV_SEP = "=";
const std::string InputParser::COMMENTS_IDEN = "[#!]";

static std::string get_last_matched(const std::string &s,
                                    const std::string &key,
                                    const std::string &vregex,
                                    int igroup)
{
    std::string sout = "";
    // a leading group to get rid of keys in comments
    std::regex r(key + SPACE_SEP + InputParser::KV_SEP + SPACE_SEP + vregex,
                 std::regex_constants::ECMAScript |
                 std::regex_constants::icase);
    std::sregex_iterator si(s.begin(), s.end(), r);
    auto ei = std::sregex_iterator();
    for (auto i = si; i != ei; i++)
    {
        std::smatch match = *i;
        // std::cout << "whole matched string: " << match.str(0) << std::endl;
        sout = match.str(igroup);
    }
    return sout;
}

void InputParser::parse_double(const std::string &vname, double &var, double de, int &flag) const
{
    flag = 0;
    std::string s = get_last_matched(params, vname,
                                     "(-?[\\d]+\\.?([\\d]+)?([ed]-?[\\d]+)?)",
                                     1);
    if (s != "")
    {
        try
        {
            var = std::stod(s);
        }
        catch (std::invalid_argument const&)
        {
            flag = 2;
        }
    }
    else
        flag = 1;
    if (flag) var = de;
}

void InputParser::parse_int(const std::string &vname, int &var, int de, int &flag) const
{
    flag = 0;
    std::string s = get_last_matched(params, vname, "(-?[\\d]+)", 1);
    if (s != "")
    {
        try
        {
            var = std::stoi(s);
        }
        catch (std::invalid_argument const&)
        {
            flag = 2;
        }
    }
    else
        flag = 1;
    if (flag) var = de;
}

void InputParser::parse_string(const std::string &vname, std::string &var, const std::string &de, int &flag) const
{
    flag = 0;
    std::string s = get_last_matched(params, vname, "([\\w\\d_\\- ,:;./]+)", 1);
    if (s != "")
        var = s;
    else
        flag = 1;
    if (flag) var = de;
}

void InputParser::parse_bool(const std::string &vname, bool &var, const bool &de, int &flag) const
{
    flag = 0;
    std::string s = get_last_matched(params, vname, "([\\w]+)", 1);
    std::transform(s.begin(), s.end(), s.begin(), ::tolower);
    if (s == "true" || s == "t" || s == ".t.")
        var = true;
    else if (s == "f" || s == "false" || s == ".f.")
        var = false;
    else
        flag = 1;
    if (flag) var = de;
}

InputParser InputFile::load(const std::string &fn, bool error_if_fail_open)
{
    std::ifstream t(fn);
    std::string params;
    if(t.is_open())
    {
        filename = fn;
        std::stringstream buffer;
        buffer << t.rdbuf();
        orig_content = buffer.str();
        // trim comment
        std::regex r(InputParser::COMMENTS_IDEN + "(.*?)\n");
        // extra \n to ensure comment in the last line is trimed
        params = std::regex_replace(orig_content + "\n", r, "\n");
    }
    else
    {
        const std::string errmsg = "Error! fail to open file " + fn;
        std::cout << errmsg << std::endl;
        if (error_if_fail_open) throw std::runtime_error(errmsg);
        std::cout << "Default parameters will be used" << std::endl;
    }
    return InputParser(params);
}

static std::string check_dirpath(const std::string &dirpath)
{
    if (dirpath.find(":") != std::string::npos)
    {
        throw std::runtime_error("input_dir contains invalid character (:) for POSIX path");
    }

    if (dirpath.back() != '/')
    {
        return dirpath + '/';
    }

    return std::string(dirpath);
}

void parse_inputfile_to_params(const std::string& fn)
{
    using namespace driver;

    InputFile inputf;
    int flag;
    auto parser = inputf.load(fn, false);

    std::string stmp;
    bool btmp;

    // TODO: invalid parameters checker

    // driver parameters
    parser.parse_string("task", driver_params.task, "rpa", flag);
    parser.parse_string("input_dir", stmp, "./", flag);
    driver_params.input_dir = check_dirpath(stmp);
    parser.parse_bool("output_gw_spec_func", driver_params.output_gw_spec_func, false, flag);

    // TODO: implement a function to read multiple double values in one line
    if (driver_params.output_gw_spec_func)
    {
        parser.parse_double("sf_omega_start", driver_params.sf_omega_start, -10.0, flag);
        parser.parse_double("sf_omega_end", driver_params.sf_omega_end, 0.0, flag);
        parser.parse_double("sf_omega_step", driver_params.sf_omega_step, 0.005, flag);
        parser.parse_int("sf_state_start", driver_params.sf_state_start, 0, flag);
        parser.parse_int("sf_state_end", driver_params.sf_state_end, 10000, flag);
        parser.parse_double("sf_gf_omega_shift", driver_params.sf_gf_omega_shift, 0.01, flag);
        parser.parse_double("sf_sigc_omega_shift", driver_params.sf_sigc_omega_shift, 0.01, flag);
    }

    // general parameters
    parser.parse_string("output_dir", stmp, "librpa.d/", flag);
    opts.set_output_dir(stmp.c_str());

    parser.parse_bool("debug", btmp, false, flag);  // backward-compatible
    if (btmp) opts.output_level = LIBRPA_VERBOSE_DEBUG;

    // chi0 related
    parser.parse_string("tfgrid_type", stmp, "minimax", flag);
    opts.tfgrids_type = get_tfgrid_type(stmp);
    parser.parse_string("parallel_routing", stmp, "auto", flag);
    std::transform(stmp.begin(), stmp.end(), stmp.begin(), ::tolower);
    opts.parallel_routing = get_parallel_routing(stmp);
    parser.parse_int("nfreq", opts.nfreq, 6, flag);

    parser.parse_bool("use_scalapack_ecrpa", btmp, false, flag);
    opts.use_scalapack_ecrpa = get_switch(btmp);

    parser.parse_bool("use_kpara_scf_eigvec", btmp, false, flag);
    opts.use_kpara_scf_eigvec = get_switch(btmp);

    parser.parse_bool("use_scalapack_gw_wc", btmp, true, flag);
    opts.use_scalapack_gw_wc = get_switch(btmp);

    parser.parse_double("cs_threshold", opts.cs_threshold, 1e-6, flag);
    parser.parse_double("vq_threshold", opts.vq_threshold, 0, flag);
    parser.parse_double("sqrt_coulomb_threshold", opts.sqrt_coulomb_threshold, 1e-4, flag);
    parser.parse_double("gf_R_threshold", opts.gf_threshold, 1e-4, flag);
    parser.parse_double("libri_chi0_threshold_C", opts.libri_chi0_threshold_C, 0.0, flag);
    parser.parse_double("libri_chi0_threshold_G", opts.libri_chi0_threshold_G, 0.0, flag);

    // exx related
    parser.parse_double("libri_exx_threshold_C", opts.libri_exx_threshold_C, 0.0, flag);
    parser.parse_double("libri_exx_threshold_D", opts.libri_exx_threshold_D, 0.0, flag);
    parser.parse_double("libri_exx_threshold_V", opts.libri_exx_threshold_V, 0.0, flag);

    // gw related
    parser.parse_double("libri_g0w0_threshold_C",  opts.libri_g0w0_threshold_C, 0.0, flag);
    parser.parse_double("libri_g0w0_threshold_G",  opts.libri_g0w0_threshold_G, 0.0, flag);
    parser.parse_double("libri_g0w0_threshold_Wc", opts.libri_g0w0_threshold_Wc, 0.0, flag);

    parser.parse_bool("replace_w_head", btmp, true, flag);
    opts.replace_w_head = get_switch(btmp);
    parser.parse_int("option_dielect_func", opts.option_dielect_func, 2, flag);

    parser.parse_bool("output_gw_sigc_mat",    btmp, false, flag);
    opts.output_gw_sigc_mat = get_switch(btmp);
    parser.parse_bool("output_gw_sigc_mat_rt", btmp, false, flag);
    opts.output_gw_sigc_mat_rt = get_switch(btmp);
    parser.parse_bool("output_gw_sigc_mat_rf", btmp, false, flag);
    opts.output_gw_sigc_mat_rf = get_switch(btmp);
}
