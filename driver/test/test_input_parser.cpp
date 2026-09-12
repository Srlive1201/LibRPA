#include "../input_parser.h"

#include <cassert>
#include <string>

void test_key_prefixes()
{
    InputParser parser(
        "tfgrids_type = minimax\n"
        "anacon_tfgrids_type = GL\n"
        "n_params_anacon_resample = 16\n"
        "input_preset = fhi-aims\n"
        "fn_basis = basis_out\n"
        "fn_basis_wfc = basis_wfc_out\n");

    int flag = 1;
    std::string value;

    parser.parse_string("tfgrids_type", value, flag);
    assert(flag == 0);
    assert(value == "minimax");

    parser.parse_string("anacon_tfgrids_type", value, flag);
    assert(flag == 0);
    assert(value == "GL");

    parser.parse_string("fn_basis", value, flag);
    assert(flag == 0);
    assert(value == "basis_out");

    parser.parse_string("input_preset", value, flag);
    assert(flag == 0);
    assert(value == "fhi-aims");

    int ivalue = 0;
    parser.parse_int("n_params_anacon_resample", ivalue, flag);
    assert(flag == 0);
    assert(ivalue == 16);
}

int main()
{
    test_key_prefixes();
}
