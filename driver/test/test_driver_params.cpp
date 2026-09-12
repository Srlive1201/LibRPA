#include "../driver.h"

#include <cassert>
#include <stdexcept>

void test_fhi_aims_preset()
{
    driver::DriverParams params;

    assert(params.input_preset == "fhi-aims");
    assert(params.fn_stru == "stru_out");
    assert(params.fn_eigocc_scf == "band_out");
    assert(params.fn_vxc_scf == "vxc_out");
    assert(params.prefix_velocity == "mommat_ks_kpt_");
}

void test_abacus_preset()
{
    driver::DriverParams params;
    params.input_preset = "abacus";
    params.apply_input_preset();

    assert(params.fn_stru == "stru_out.txt");
    assert(params.fn_eigocc_scf == "band_out.txt");
    assert(params.fn_vxc_scf == "vxc_out.txt");
    assert(params.prefix_velocity == "velocity_matrix");
    assert(params.fn_basis_wfc == "basis_wfc_out");
}

void test_invalid_preset()
{
    driver::DriverParams params;
    params.input_preset = "unknown";

    bool threw = false;
    try
    {
        params.apply_input_preset();
    }
    catch (const std::invalid_argument &)
    {
        threw = true;
    }
    assert(threw);
}

int main()
{
    test_fhi_aims_preset();
    test_abacus_preset();
    test_invalid_preset();
}
