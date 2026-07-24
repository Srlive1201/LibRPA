#ifdef NDEBUG
#undef NDEBUG
#endif

#include "../driver.h"
#include "../inputfile.h"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <stdexcept>
#include <string>

namespace
{

class TemporaryInput
{
public:
    explicit TemporaryInput(const std::string& contents)
        : path_("test_qsgw_inputfile.tmp")
    {
        std::ofstream output(path_, std::ios::trunc);
        output << contents;
    }

    ~TemporaryInput() { std::remove(path_.c_str()); }
    const std::string& path() const { return path_; }

private:
    std::string path_;
};

void reset()
{
    driver::driver_params = driver::DriverParams{};
    driver::opts = librpa::Options{};
    driver::n_spinor = 1;
}

void parse(const std::string& contents)
{
    reset();
    const TemporaryInput input(contents);
    parse_inputfile_to_params(input.path());
}

template <typename Function>
void assert_throws(Function&& function)
{
    bool threw = false;
    try
    {
        function();
    }
    catch (const std::exception&)
    {
        threw = true;
    }
    assert(threw);
}

template <typename Function>
void assert_throws_with_message(Function&& function,
                                const std::string& expected)
{
    bool threw = false;
    try
    {
        function();
    }
    catch (const std::exception& error)
    {
        threw = true;
        assert(std::string(error.what()).find(expected) !=
               std::string::npos);
    }
    assert(threw);
}

std::string valid_qsgw_prefix()
{
    return "task = qsgw\n";
}

std::string valid_qsgw_band_prefix()
{
    return "task = qsgw_band\n";
}

void test_defaults_and_explicit_linear_mixing()
{
    const driver::DriverParams defaults;
    assert(defaults.qsgw_input_contract == "qsgw_input.contract");
    assert(defaults.qsgw_mixer == "none");
    assert(std::abs(defaults.qsgw_mixing_beta - 0.2) < 1.0e-14);
    assert(defaults.qsgw_min_iter == 1);
    assert(defaults.qsgw_max_iter == 10);
    assert(defaults.qsgw_band0_unoccupied_keep == 10);
    assert(defaults.qsgw_band0_cut_mode == 0);
    assert(std::abs(defaults.qsgw_band0_cut_shift_ha - 20.0) < 1.0e-14);

    parse(valid_qsgw_prefix() +
          "qsgw_input_contract = manifests/si.contract\n"
          "qsgw_mixer = NONE\n"
          "qsgw_mixing_beta = 3.5D-1\n"
          "qsgw_min_iter = 5\n"
          "qsgw_max_iter = 10\n"
          "qsgw_write_iteration_matrices = true\n"
          "qsgw_convergence_tolerance_ev = 2e-5\n");
    assert(driver::driver_params.qsgw_input_contract ==
           "manifests/si.contract");
    assert(driver::driver_params.qsgw_mixer == "none");
    assert(std::abs(driver::driver_params.qsgw_mixing_beta - 0.35) <
           1.0e-14);
    assert(driver::driver_params.qsgw_min_iter == 5);
    assert(driver::driver_params.qsgw_max_iter == 10);
    assert(driver::driver_params.qsgw_write_iteration_matrices);
    const std::string formatted = driver::driver_params.format();
    assert(formatted.find("qsgw_input_contract") != std::string::npos);
    assert(formatted.find("qsgw_write_iteration_matrices = true") !=
           std::string::npos);
}

void test_qsgw_accepts_only_same_grid_head_only()
{
    for (const std::string& prefix :
         {valid_qsgw_prefix(), valid_qsgw_band_prefix()})
    {
        parse(prefix +
              "replace_w_head = true\n"
              "option_dielect_func = 4\n");
        assert(driver::opts.replace_w_head == LIBRPA_SWITCH_ON);
        assert(driver::opts.option_dielect_func == 4);
        assert(!driver::driver_params.use_pyatb);

        assert_throws_with_message([&] {
            parse(prefix +
                  "replace_w_head = true\n"
                  "option_dielect_func = 3\n");
        }, "QSGW supports only analytic head-only mode");
        assert_throws_with_message([&] {
            parse(prefix +
                  "replace_w_head = true\n"
                  "option_dielect_func = 4\n"
                  "use_pyatb = true\n");
        }, "QSGW independent PyATB head updates are unsupported");
        assert_throws_with_message([&] {
            parse(prefix + "use_pyatb = true\n");
        }, "QSGW independent PyATB head updates are unsupported");
    }
}

void test_qsgw_requires_replicated_scf_wavefunctions()
{
    const std::string message =
        "QSGW fixed-basis iteration requires replicated SCF wavefunctions";
    for (const std::string& prefix :
         {valid_qsgw_prefix(), valid_qsgw_band_prefix()})
    {
        assert_throws_with_message([&] {
            parse(prefix + "use_kpara_scf_eigvec = true\n");
        }, message);
    }
}

void test_qsgw_band_uses_the_qsgw_contract_and_iteration_controls()
{
    parse(valid_qsgw_band_prefix() +
          "qsgw_input_contract = manifests/si-band.contract\n"
          "qsgw_mixer = LINEAR\n"
          "qsgw_mixing_beta = 0.27\n"
          "qsgw_min_iter = 2\n"
          "qsgw_max_iter = 5\n"
          "qsgw_band0_unoccupied_keep = 7\n"
          "qsgw_band0_cut_mode = 1\n"
          "qsgw_band0_cut_shift_ha = 3.5\n"
          "qsgw_write_iteration_matrices = true\n");
    assert(driver::driver_params.qsgw_input_contract ==
           "manifests/si-band.contract");
    assert(driver::driver_params.qsgw_mixer == "linear");
    assert(std::abs(driver::driver_params.qsgw_mixing_beta - 0.27) <
           1.0e-14);
    assert(driver::driver_params.qsgw_min_iter == 2);
    assert(driver::driver_params.qsgw_max_iter == 5);
    assert(driver::driver_params.qsgw_band0_unoccupied_keep == 7);
    assert(driver::driver_params.qsgw_band0_cut_mode == 1);
    assert(std::abs(driver::driver_params.qsgw_band0_cut_shift_ha - 3.5) <
           1.0e-14);
    assert(driver::driver_params.qsgw_write_iteration_matrices);
    assert(driver::driver_params.format().find(
               "qsgw_input_contract = manifests/si-band.contract") !=
           std::string::npos);
    assert(driver::driver_params.format().find(
               "qsgw_band0_cut_mode = 1") != std::string::npos);
}

void test_qsgw_does_not_require_the_g0w0_matrix_dump_switch()
{
    parse("task = qsgw\noutput_gw_sigc_ks_mat_kf = false\n");
    assert(driver::opts.output_gw_sigc_ks_mat_kf == LIBRPA_SWITCH_OFF);

    parse("task = qsgw_band\n");
    assert(driver::opts.output_gw_sigc_ks_mat_kf == LIBRPA_SWITCH_OFF);
}

void test_qsgw_preserves_upstream_crystal_symmetry_options()
{
    parse(valid_qsgw_prefix() + "use_symmetry_exx = true\n");
    assert(driver::opts.use_symmetry_exx == LIBRPA_SWITCH_ON);
    assert(driver::opts.use_symmetry_gw == LIBRPA_SWITCH_OFF);
    assert(driver::opts.use_symmetry_rpa == LIBRPA_SWITCH_OFF);

    parse(valid_qsgw_prefix() + "use_symmetry_gw = true\n");
    assert(driver::opts.use_symmetry_exx == LIBRPA_SWITCH_OFF);
    assert(driver::opts.use_symmetry_gw == LIBRPA_SWITCH_ON);
    assert(driver::opts.use_symmetry_rpa == LIBRPA_SWITCH_ON);

    parse(valid_qsgw_prefix() + "use_symmetry_rpa = true\n");
    assert(driver::opts.use_symmetry_exx == LIBRPA_SWITCH_OFF);
    assert(driver::opts.use_symmetry_gw == LIBRPA_SWITCH_OFF);
    assert(driver::opts.use_symmetry_rpa == LIBRPA_SWITCH_ON);

    parse(valid_qsgw_prefix() +
          "use_symmetry_exx = true\n"
          "use_symmetry_gw = true\n"
          "use_symmetry_rpa = true\n");
    assert(driver::opts.use_symmetry_exx == LIBRPA_SWITCH_ON);
    assert(driver::opts.use_symmetry_gw == LIBRPA_SWITCH_ON);
    assert(driver::opts.use_symmetry_rpa == LIBRPA_SWITCH_ON);
}

void test_staged_qsgw_rejects_ambiguous_or_unsupported_inputs()
{
    for (const std::string& removed : {
             "qsgw_update_hartree = false\n",
             "qsgw_hartree_coulomb = full\n",
             "qsgw_hartree_normalization = weighted_occupations\n",
             "qsgw_export_hamiltonian_for_pyatb = false\n",
             "qsgw_hr_export_full_mp_rgrid = false\n"})
    {
        assert_throws_with_message([&] {
            parse(valid_qsgw_prefix() + removed);
        }, "is not supported by the head-only QSGW workflow");
    }
    assert_throws([&] {
        parse(valid_qsgw_prefix() + "qsgw_mixer = pulay\n");
    });
    assert_throws([&] {
        parse(valid_qsgw_prefix() + "qsgw_mixing_beta = 0\n");
    });
    assert_throws([&] {
        parse(valid_qsgw_prefix() +
              "qsgw_min_iter = 5\nqsgw_max_iter = 4\n");
    });
    assert_throws([&] {
        parse(valid_qsgw_band_prefix() +
              "qsgw_band0_unoccupied_keep = -1\n");
    });
    assert_throws([&] {
        parse(valid_qsgw_band_prefix() +
              "qsgw_band0_cut_mode = 3\n");
    });
    assert_throws([&] {
        parse(valid_qsgw_band_prefix() +
              "qsgw_band0_cut_shift_ha = nan\n");
    });
    assert_throws_with_message([&] {
        parse(valid_qsgw_prefix() +
              "qsgw_max_iter = ten\n");
    }, "qsgw_max_iter must be a valid integer");
    assert_throws_with_message([&] {
        parse(valid_qsgw_prefix() +
              "qsgw_mixing_beta = 0.2garbage\n");
    }, "qsgw_mixing_beta must be a valid floating-point value");
    assert_throws_with_message([&] {
        parse(valid_qsgw_prefix() +
              "qsgw_max_iter = 5\n"
              "qsgw_max_iter = invalid\n");
    }, "qsgw_max_iter must be a valid integer");
    assert_throws([&] {
        parse(valid_qsgw_prefix() +
              "replace_w_head = true\noption_dielect_func = 0\n");
    });
    assert_throws([&] {
        parse(valid_qsgw_prefix() +
              "replace_w_head = true\noption_dielect_func = 2\n");
    });
}

void test_g0w0_parser_defaults_are_unchanged()
{
    parse(
        "task = g0w0\n"
        "qsgw_mixer = invalid\n"
        "qsgw_mixing_beta = -1\n"
        "qsgw_band0_unoccupied_keep = -1\n"
        "qsgw_band0_cut_mode = 99\n"
        "qsgw_band0_cut_shift_ha = nan\n"
        "use_symmetry_exx = true\n"
        "use_symmetry_gw = true\n"
        "use_symmetry_rpa = true\n");
    assert(driver::driver_params.qsgw_mixer == "none");
    assert(std::abs(driver::driver_params.qsgw_mixing_beta - 0.2) <
           1.0e-14);
    assert(driver::driver_params.qsgw_band0_unoccupied_keep == 10);
    assert(driver::driver_params.qsgw_band0_cut_mode == 0);
    assert(std::abs(driver::driver_params.qsgw_band0_cut_shift_ha - 20.0) <
           1.0e-14);
    assert(driver::opts.output_gw_sigc_ks_mat_kf == LIBRPA_SWITCH_OFF);
    assert(driver::opts.use_symmetry_exx == LIBRPA_SWITCH_ON);
    assert(driver::opts.use_symmetry_gw == LIBRPA_SWITCH_ON);
    assert(driver::opts.use_symmetry_rpa == LIBRPA_SWITCH_ON);
    assert(driver::driver_params.format().find("qsgw_") ==
           std::string::npos);

    parse(
        "task = g0w0\n"
        "replace_w_head = true\n"
        "option_dielect_func = 3\n"
        "use_pyatb = true\n");
    assert(driver::opts.replace_w_head == LIBRPA_SWITCH_ON);
    assert(driver::opts.option_dielect_func == 3);
    assert(driver::driver_params.use_pyatb);
}

} // namespace

int main()
{
    test_defaults_and_explicit_linear_mixing();
    test_qsgw_accepts_only_same_grid_head_only();
    test_qsgw_requires_replicated_scf_wavefunctions();
    test_qsgw_band_uses_the_qsgw_contract_and_iteration_controls();
    test_qsgw_does_not_require_the_g0w0_matrix_dump_switch();
    test_qsgw_preserves_upstream_crystal_symmetry_options();
    test_staged_qsgw_rejects_ambiguous_or_unsupported_inputs();
    test_g0w0_parser_defaults_are_unchanged();
    return 0;
}
