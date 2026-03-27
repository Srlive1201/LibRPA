#include "inputfile.h"

#include <regex>
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>

#include "driver.h"

#include "librpa_enums.h"
#include "../src/utils/constants.h"

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

void InputParser::parse_double(const std::string &vname, double &var, int &flag) const
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
}


void InputParser::parse_double(const std::string &vname, double &var, double de, int &flag) const
{
    parse_double(vname, var, flag);
    if (flag) var = de;
}

void InputParser::parse_int(const std::string &vname, int &var, int &flag) const
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
}

void InputParser::parse_int(const std::string &vname, int &var, int de, int &flag) const
{
    parse_int(vname, var, flag);
    if (flag) var = de;
}

void InputParser::parse_string(const std::string &vname, std::string &var, int &flag) const
{
    flag = 0;
    std::string s = get_last_matched(params, vname, "([\\w\\d_\\- ,:;./]+)", 1);
    if (s != "")
        var = s;
    else
        flag = 1;
}

void InputParser::parse_string(const std::string &vname, std::string &var, const std::string &de, int &flag) const
{
    parse_string(vname, var, flag);
    if (flag) var = de;
}

void InputParser::parse_bool(const std::string &vname, bool &var, int &flag) const
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
}

void InputParser::parse_bool(const std::string &vname, bool &var, const bool &de, int &flag) const
{
    parse_bool(vname, var, flag);
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

#define _parse_int(obj, name) parser.parse_int(#name, obj.name, flag)
#define _parse_double(obj, name) parser.parse_double(#name, obj.name, flag)
#define _parse_switch(obj, name) parser.parse_bool(#name, btmp, flag); if (flag == 0) obj.name = get_switch(btmp);
#define _parse_string_post(obj, name, post) parser.parse_string(#name, stmp, flag); if (flag == 0) obj.name = post(stmp);

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
    parser.parse_string("task", driver_params.task, "unset", flag);
    parser.parse_string("constants_choice", driver_params.constants_choice, "internal", flag);
    if (driver_params.constants_choice == "aims")
    {
        librpa_int::set_aims_constants();
    }
    _parse_string_post(driver_params, input_dir, check_dirpath);
    _parse_double(driver_params, cs_threshold);
    parser.parse_bool("output_gw_spec_func", driver_params.output_gw_spec_func, false, flag);

    // TODO: implement a function to read multiple double values in one line
    if (driver_params.output_gw_spec_func)
    {
        _parse_double(driver_params, sf_omega_start);
        _parse_double(driver_params, sf_omega_end);
        _parse_double(driver_params, sf_omega_step);
        _parse_int(driver_params, sf_state_start);
        _parse_int(driver_params, sf_state_end);
        _parse_double(driver_params, sf_gf_omega_shift);
        _parse_double(driver_params, sf_sigc_omega_shift);
    }

    // general parameters
    parser.parse_string("output_dir", stmp, "librpa.d/", flag);
    opts.set_output_dir(stmp.c_str());
    _parse_string_post(opts, parallel_routing, get_parallel_routing);

    parser.parse_bool("debug", btmp, false, flag);  // backward-compatible
    if (btmp) opts.output_level = LIBRPA_VERBOSE_DEBUG;
    parser.parse_string("output_level", stmp, "info", flag);
    if (flag == 0) opts.output_level = get_verbose(stmp);
    _parse_double(opts, vq_threshold);
    _parse_switch(opts, use_kpara_scf_eigvec);

    _parse_string_post(opts, tfgrids_type, get_tfgrid_type);
    if (flag != 0) // backward compatible
    {
        parser.parse_string("tfgrid_type", stmp, flag);
        if (flag == 0) opts.tfgrids_type = get_tfgrid_type(stmp);
    }
    if (opts.tfgrids_type == LibrpaTimeFreqGrid::TFGRID_UNSET)
        opts.tfgrids_type = LibrpaTimeFreqGrid::Minimax;
    _parse_int(opts, nfreq);

    // RPA specific
    parser.parse_double("gf_R_threshold", opts.gf_threshold, flag); // backward compatible
    _parse_double(opts, gf_threshold);
    _parse_double(opts, libri_chi0_threshold_C);
    _parse_double(opts, libri_chi0_threshold_G);
    _parse_switch(opts, use_scalapack_ecrpa);

    // EXX specific
    _parse_double(opts, libri_exx_threshold_C);
    _parse_double(opts, libri_exx_threshold_D);
    _parse_double(opts, libri_exx_threshold_V);

    // GW specific
    _parse_double(opts, sqrt_coulomb_threshold);
    _parse_switch(opts, use_scalapack_gw_wc);
    _parse_double(opts, libri_g0w0_threshold_C);
    _parse_double(opts, libri_g0w0_threshold_G);
    _parse_double(opts, libri_g0w0_threshold_Wc);
    _parse_switch(opts, replace_w_head);
    _parse_int(opts, option_dielect_func);
    _parse_switch(opts, output_gw_sigc_mat);
    _parse_switch(opts, output_gw_sigc_mat_rt);
    _parse_switch(opts, output_gw_sigc_mat_rf);
}

#undef _parse_int
#undef _parse_double
#undef _parse_switch
#undef _parse_string_post
