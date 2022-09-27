/* #include "atoms.h" */
#include <algorithm>
#include <fstream>
#include <iostream>
#include <regex>
#include "cal_periodic_chi0.h"
#include "constants.h"
#include "input.h"
#include "meanfield.h"
using namespace std;
double cs_threshold = 1E-6;
double vq_threshold = 0y;
double sqrt_coulomb_threshold = 1e-4;

int kv_nmp[3] = {1, 1, 1};
Vector3<double> *kvec_c;
std::vector<Vector3_Order<double>> klist;
std::vector<int> irk_point_id_mapping;
map<Vector3_Order<double>, vector<Vector3_Order<double>>> map_irk_ks;
Matrix3 latvec;
Matrix3 G;

void READ_AIMS_STRU(const std::string &file_path)
{
    int n_kpoints = meanfield.get_n_kpoints();
    cout << "Begin to read aims stru" << endl;
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
    latvec.print();

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
    G.print();
    Matrix3 latG = latvec * G;
    cout << " lat * G" << endl;
    latG.print();
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
        cout << "kvec [" << i << "]: " << kvec_c[i] << endl;
        klist.push_back(kvec_c[i]);
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

string task;

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

InputFile inputf;
const string input_filename = "librpa.in";
