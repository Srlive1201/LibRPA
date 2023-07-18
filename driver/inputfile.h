#pragma once
#include <string>
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

void parse_inputfile_to_params(const std::string& fn);

extern const std::string input_filename;
