/* #include "atoms.h" */
#include <algorithm>
#include <fstream>
#include <iostream>
#include <regex>
#include "constants.h"
#include "input.h"
using namespace std;

int kv_nmp[3] = {1, 1, 1};
Vector3<double> *kvec_c;
std::vector<Vector3_Order<double>> klist;
std::vector<Vector3_Order<double>> kfrac_list;
std::vector<int> irk_point_id_mapping;
map<Vector3_Order<double>, vector<Vector3_Order<double>>> map_irk_ks;
Matrix3 latvec;
std::array<std::array<double, 3>, 3> lat_array;
Matrix3 G;

void READ_AIMS_STRU(const int& n_kpoints, const std::string &file_path)
{
    // cout << "Begin to read aims stru" << endl;
    ifstream infile;
    string x, y, z;
    infile.open(file_path);
    infile >> x >> y >> z;
    latvec.e11 = stod(x);
    latvec.e12 = stod(y);
    latvec.e13 = stod(z);
    infile >> x >> y >> z;
    latvec.e21 = stod(x);
    latvec.e22 = stod(y);
    latvec.e23 = stod(z);
    infile >> x >> y >> z;
    latvec.e31 = stod(x);
    latvec.e32 = stod(y);
    latvec.e33 = stod(z);
    latvec /= ANG2BOHR;
    // latvec.print();
    lat_array[0] = {latvec.e11,latvec.e12,latvec.e13};
    lat_array[1] = {latvec.e21,latvec.e22,latvec.e23};
    lat_array[2] = {latvec.e31,latvec.e32,latvec.e33};

    infile >> x >> y >> z;
    G.e11 = stod(x);
    G.e12 = stod(y);
    G.e13 = stod(z);
    infile >> x >> y >> z;
    G.e21 = stod(x);
    G.e22 = stod(y);
    G.e23 = stod(z);
    infile >> x >> y >> z;
    G.e31 = stod(x);
    G.e32 = stod(y);
    G.e33 = stod(z);

    G /= TWO_PI;
    G *= ANG2BOHR;
    // G.print();
    // Matrix3 latG = latvec * G.Transpose();
    // cout << " lat * G^T" << endl;
    // latG.print();
    infile >> x >> y >> z;
    kv_nmp[0] = stoi(x);
    kv_nmp[1] = stoi(y);
    kv_nmp[2] = stoi(z);
    kvec_c = new Vector3<double>[n_kpoints];
    for (int i = 0; i != n_kpoints; i++)
    {
        infile >> x >> y >> z;
        kvec_c[i] = {stod(x), stod(y), stod(z)};
        kvec_c[i] *= (ANG2BOHR / TWO_PI);
        // cout << "kvec [" << i << "]: " << kvec_c[i] << endl;
        Vector3_Order<double> k(kvec_c[i]);
        klist.push_back(k);
        kfrac_list.push_back(latvec * k);
    }
    for (int i = 0; i != n_kpoints; i++)
    {
        infile >> x;
        int id_irk = stoi(x) - 1;
        irk_point_id_mapping.push_back(id_irk);
        map_irk_ks[klist[id_irk]].push_back(klist[i]);
    }
}

const std::string SPACE_SEP = "[ \r\f\t]*";
const std::string InputParser::KV_SEP = "=";
const std::string InputParser::COMMENTS_IDEN = "[#!]";

std::string get_last_matched(const std::string &s,
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
    std::string s = get_last_matched(params, vname, "([\\w ,;.]+)", 1);
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
        const string errmsg = "Error! fail to open file " + fn;
        std::cout << errmsg << std::endl;
        if (error_if_fail_open) throw runtime_error(errmsg);
        std::cout << "Default parameters will be used" << std::endl;
    }
    return InputParser(params);
}

void parse_inputfile_to_params(const string& fn)
{
    InputFile inputf;
    int flag;
    auto parser = inputf.load(fn, false);

    // general parameters
    parser.parse_string("task", Params::task, "rpa", flag);
    parser.parse_bool("debug", Params::debug, false, flag);

    // chi0 related
    parser.parse_string("tfgrid_type", Params::tfgrids_type, "minimax", flag);
    parser.parse_string("chi_parallel_routing", Params::chi_parallel_routing, "auto", flag);
    parser.parse_int("nfreq", Params::nfreq, 6, flag);
    parser.parse_bool("use_libri_chi0", Params::use_libri_chi0, false, flag);
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
    parser.parse_bool("use_libri_exx", Params::use_libri_exx, false, flag);
    parser.parse_string("exx_parallel_routing", Params::exx_parallel_routing, "auto", flag);
    parser.parse_double("libri_exx_threshold_CSM", Params::libri_exx_threshold_CSM, 0.0, flag);
    parser.parse_double("libri_exx_threshold_C", Params::libri_exx_threshold_C, 0.0, flag);
    parser.parse_double("libri_exx_threshold_D", Params::libri_exx_threshold_D, 0.0, flag);
    parser.parse_double("libri_exx_threshold_V", Params::libri_exx_threshold_V, 0.0, flag);

    // gw related
    parser.parse_bool("use_libri_gw", Params::use_libri_gw, true, flag);
    parser.parse_bool("output_gw_sigc_mat", Params::output_gw_sigc_mat, true, flag);
}

const string input_filename = "librpa.in";
