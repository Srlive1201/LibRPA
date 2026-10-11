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

    params.input_preset = "abacus";
    params.apply_input_preset();
    params.input_preset = "aims";
    params.apply_input_preset();

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
    assert(params.fn_bz_sampling == "bz_sample.txt");
    assert(params.fn_basis == "basis_map.txt");
    assert(params.fn_basis_wfc == "wfc_basis.txt");
    assert(params.fn_basis_aux == "aux_basis.txt");
    assert(params.fn_basis_aux_shrink == "aux_basis_s.txt");
    assert(params.fn_eigocc_scf == "band_out.txt");
    assert(params.fn_vxc_scf == "vxc.txt");
    assert(params.prefix_lri_coeff == "Cs_");
    assert(params.prefix_lri_coeff_shrink == "Cs_shrink_");
    assert(params.prefix_shrink_sinvS == "sinvS_");
    assert(params.prefix_coul_full == "V_full_");
    assert(params.prefix_coul_cut == "V_cut_");
    assert(params.prefix_eigvecs_scf == "KS_wfc_");
    assert(params.prefix_velocity == "velocity_matrix");

    params.input_preset = "abacus-legacy";
    params.apply_input_preset();

    assert(params.input_preset == "abacus-legacy");
    assert(params.fn_stru == "stru_out");
    assert(params.fn_eigocc_scf == "band_out");
    assert(params.fn_vxc_scf == "vxc_out");
    assert(params.prefix_velocity == "velocity_matrix");
    assert(params.fn_basis_wfc == "basis_wfc_out");
    assert(params.prefix_lri_coeff == "Cs_data");
    assert(params.prefix_lri_coeff_shrink == "Cs_shrinked_data");
    assert(params.prefix_shrink_sinvS == "shrink_sinvS_");
    assert(params.prefix_coul_full == "coulomb_mat");
    assert(params.prefix_coul_cut == "coulomb_cut");
    assert(params.prefix_eigvecs_scf == "KS_eigenvector");
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
