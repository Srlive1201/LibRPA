#include "inputfile.h"

#include <regex>
#include <fstream>
#include <sstream>
#include <iostream>

#include "params.h"

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

void InputParser::parse_double(const std::string &vname, double &var, double de, int &flag)
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
        catch (std::invalid_argument)
        {
            flag = 2;
        }
    }
    else
        flag = 1;
    if (flag) var = de;
}

void InputParser::parse_int(const std::string &vname, int &var, int de, int &flag)
{
    flag = 0;
    std::string s = get_last_matched(params, vname, "(-?[\\d]+)", 1);
    if (s != "")
    {
        try
        {
            var = std::stoi(s);
        }
        catch (std::invalid_argument)
        {
            flag = 2;
        }
    }
    else
        flag = 1;
    if (flag) var = de;
}

void InputParser::parse_string(const std::string &vname, std::string &var, const std::string &de, int &flag)
{
    flag = 0;
    std::string s = get_last_matched(params, vname, "([\\w\\d_\\- ,;.]+)", 1);
    if (s != "")
        var = s;
    else
        flag = 1;
    if (flag) var = de;
}

void InputParser::parse_bool(const std::string &vname, bool &var, const bool &de, int &flag)
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

void parse_inputfile_to_params(const std::string& fn)
{
    InputFile inputf;
    int flag;
    auto parser = inputf.load(fn, false);

    // general parameters
    parser.parse_string("task", Params::task, "rpa", flag);
    parser.parse_string("output_dir", Params::output_dir, "librpa.d", flag);
    parser.parse_string("output_file", Params::output_file, "stdout", flag);
    parser.parse_bool("debug", Params::debug, false, flag);
    parser.parse_bool("binary_input", Params::binary_input, false, flag);

    // chi0 related
    parser.parse_string("tfgrid_type", Params::tfgrids_type, "minimax", flag);
    parser.parse_string("parallel_routing", Params::parallel_routing, "auto", flag);
    parser.parse_int("nfreq", Params::nfreq, 6, flag);
    parser.parse_bool("use_scalapack_ecrpa", Params::use_scalapack_ecrpa, false, flag);
    parser.parse_bool("use_scalapack_gw_wc", Params::use_scalapack_gw_wc, false, flag);
    parser.parse_double("cs_threshold", Params::cs_threshold, 1e-6, flag);
    parser.parse_double("vq_threshold", Params::vq_threshold, 0, flag);
    parser.parse_double("sqrt_coulomb_threshold", Params::sqrt_coulomb_threshold, 1e-4, flag);
    parser.parse_double("gf_R_threshold", Params::gf_R_threshold, 1e-4, flag);
    parser.parse_double("libri_chi0_threshold_CSM", Params::libri_chi0_threshold_CSM, 0.0, flag);
    parser.parse_double("libri_chi0_threshold_C", Params::libri_chi0_threshold_C, 0.0, flag);
    parser.parse_double("libri_chi0_threshold_G", Params::libri_chi0_threshold_G, 0.0, flag);

    // exx related
    parser.parse_double("libri_exx_threshold_CSM", Params::libri_exx_threshold_CSM, 0.0, flag);
    parser.parse_double("libri_exx_threshold_C", Params::libri_exx_threshold_C, 0.0, flag);
    parser.parse_double("libri_exx_threshold_D", Params::libri_exx_threshold_D, 0.0, flag);
    parser.parse_double("libri_exx_threshold_V", Params::libri_exx_threshold_V, 0.0, flag);

    // gw related
    parser.parse_bool("output_gw_sigc_mat", Params::output_gw_sigc_mat, true, flag);
    parser.parse_bool("replace_w_head", Params::replace_w_head, true, flag);
    parser.parse_int("option_dielect_func", Params::option_dielect_func, 2, flag);
}

const std::string input_filename = "librpa.in";
