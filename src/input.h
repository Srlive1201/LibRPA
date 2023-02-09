/*
 * @file input.h
 */
#ifndef INPUT_H
#define INPUT_H

/* #include "atoms.h" */
#include <array>
#include <string>
#include <vector>
#include <map>
#include "params.h"
#include "vector3.h"
#include "vector3_order.h"
#include "matrix3.h"

extern int kv_nmp[3];
//! lattice vectors as a 3D-matrix, each row as a lattice vector
extern Matrix3 latvec;
//! same as latvec, but a nested array for LibRI call
extern std::array<std::array<double, 3>, 3> lat_array;
//! reciprocal lattice vectors as a 3D-matrix, each row as a reciprocal vector
extern Matrix3 G;
extern std::vector<Vector3_Order<double>> klist;
extern std::vector<Vector3_Order<double>> kfrac_list;
extern std::vector<int> irk_point_id_mapping;
extern map<Vector3_Order<double>, vector<Vector3_Order<double>>> map_irk_ks;
extern Vector3<double> *kvec_c;

void READ_AIMS_STRU(const int& n_kpoints, const std::string &file_path);

//! Class to parse parameters from string
/*!
  Member functions prefixed with parse_ have a int argument at last.
  This flag is to indicate whether parsing is succesful.
  The meaning of return value:

  - 0: succesful
  - 1: failed to find the parameter string, or the parameter value is empty
  - 2: find the parameter string, but fail to convert the value to requested type

  When multiple matches are found, the last one will be used.

  Only successful parse will lead to value change of the variable.
  Otherwise the default value will be used.
 */
class InputParser
{
    public:
        //! separator between the parameter key and value
        static const std::string KV_SEP;
        //! comments identifier to trim
        static const std::string COMMENTS_IDEN;
    private:
        std::string params; //!< string storing all parameters, comments trimed
    public:
        InputParser() {};
        InputParser(const std::string &ps): params(ps) {};
        std::string get_params() const { return params; };
        void load(const std::string &ps) { params = ps; };
        //! parse single double value to var
        void parse_double(const std::string &vname, double &var, double de, int &flag);
        //! parse single integer value to var
        void parse_int(const std::string &vname, int &var, int de, int &flag);
        //! parse single string to var
        void parse_string(const std::string &vname, std::string &var, const std::string &de, int &flag);
        //! parse single bool to var
        void parse_bool(const std::string &vname, bool &var, const bool &de, int &flag);
};

//! Class to represent the input file
/*!
  Use the member function load to get the parameter parser.
 */
class InputFile
{
    private:
        std::string filename;
        std::string orig_content;
    public:
        InputFile() {};
        InputParser load(const std::string &fn, bool error_if_fail_open = false);
        std::string get_filename() { return filename; };
        std::string get_orig_content() { return orig_content; };
};

void parse_inputfile_to_params(const string& fn, Params& p);

extern const string input_filename;
#endif
