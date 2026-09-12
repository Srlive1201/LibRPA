#include "inputfile.h"

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <cstring>
#include <regex>
#include <sstream>
#include <limits>

#include "driver.h"

#include <librpa_enums.h>

#include "../src/io/fs.h"
#include "../src/utils/constants.h"
#include "../src/utils/dev_options.h"

InputParser InputFile::load(const std::string &fn, bool error_if_fail_open)
{
    if (error_if_fail_open)
    {
        librpa_int::require_readable_file(fn);
    }
    std::ifstream t(fn);
    std::string params;
    if (t.is_open())
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
#define _parse_bool(obj, name) parser.parse_bool(#name, obj.name, flag)
#define _parse_string(obj, name) parser.parse_string(#name, obj.name, flag);
#define _parse_switch(obj, name) parser.parse_bool(#name, btmp, flag); if (flag == 0) obj.name = get_switch(btmp);
#define _parse_string_post(obj, name, post) parser.parse_string(#name, stmp, flag); if (flag == 0) obj.name = post(stmp);

static void validate_input_parameters()
{
    const auto &params = driver::driver_params;
    if (params.output_gw_spec_func)
    {
        if (params.sf_omega_step <= 0.0)
            throw std::runtime_error("sf_omega_step must be positive");
        if (params.sf_omega_end < params.sf_omega_start)
            throw std::runtime_error("sf_omega_end must be no smaller than sf_omega_start");
        if (params.sf_state_start >= 0 && params.sf_state_end >= 0
            && params.sf_state_end <= params.sf_state_start)
            throw std::runtime_error("sf_state_end must be greater than sf_state_start");
    }
}

void parse_inputfile_to_params(const std::string &fn)
{
    using namespace driver;
    using librpa_int::global::dev_opts;

    InputFile inputf;
    int flag;
    auto parser = inputf.load(fn, true);

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
    parser.parse_string("input_preset", driver_params.input_preset, "fhi-aims", flag);
    driver_params.apply_input_preset();
    _parse_string_post(driver_params, input_dir, check_dirpath);
    _parse_double(driver_params, cs_threshold);
    _parse_bool(driver_params, output_energy_qp);
    _parse_int(driver_params, i_state_low);
    _parse_int(driver_params, i_state_high);
    _parse_bool(driver_params, output_hamgnn);
    _parse_bool(driver_params, use_pyatb);
    _parse_bool(driver_params, output_gw_spec_func);
    _parse_bool(driver_params, use_spinor_wfc);
    parser.parse_bool("use_soc", driver_params.use_spinor_wfc, flag);  // backward-compatible
    if (driver_params.use_spinor_wfc)
    {
        driver::n_spinor = 2;
    }
    _parse_string(driver_params, prefix_lri_coeff);
    _parse_string(driver_params, prefix_lri_coeff_shrink);
    _parse_string(driver_params, prefix_shrink_sinvS);
    _parse_string(driver_params, prefix_coul_full);
    _parse_string(driver_params, prefix_coul_cut);
    _parse_string(driver_params, prefix_eigvecs_scf);
    _parse_string(driver_params, fn_stru);
    _parse_string(driver_params, fn_basis);
    _parse_string(driver_params, fn_basis_wfc);
    _parse_string(driver_params, fn_basis_aux);
    _parse_string(driver_params, fn_basis_aux_shrink);
    parser.parse_string("fn_basis_shrink", driver_params.fn_basis_aux_shrink, flag);
    _parse_string(driver_params, fn_bz_sampling);
    _parse_string(driver_params, fn_eigocc_scf);
    _parse_string(driver_params, fn_dielfunc);
    _parse_string(driver_params, fn_vxc_scf);
    _parse_string(driver_params, prefix_velocity);
    _parse_string(driver_params, fn_band_kpath_info);
    _parse_int(driver_params, version_coul_reader);
    _parse_int(driver_params, version_lri_reader);

    // TODO: implement a function to read multiple double values in one line
    if (driver_params.output_gw_spec_func)
    {
        _parse_double(driver_params, sf_omega_start);
        _parse_double(driver_params, sf_omega_end);
        _parse_double(driver_params, sf_omega_step);
        _parse_int(driver_params, sf_state_start);
        _parse_int(driver_params, sf_state_end);
    }

    // Development options, see src/utils/dev_options.h
    // Only for prototyping, testing and debugging in C++
    _parse_bool(dev_opts, use_chi0_q_uhap_split);
    _parse_bool(dev_opts, use_delayed_ft_shrink);
    _parse_bool(dev_opts, use_real_dense_gw_wc);

    // general runtime parameters
    parser.parse_string("output_dir", stmp, flag);
    if (flag == 0)
        opts.set_output_dir(stmp.c_str());
    parser.parse_string("restart_from_dir", stmp, flag);
    if (flag == 0)
        opts.set_restart_from_dir(stmp.c_str());
    _parse_string_post(opts, parallel_routing, get_parallel_routing);

    parser.parse_bool("debug", btmp, false, flag);  // backward-compatible
    if (btmp) driver_params.output_level = LIBRPA_VERBOSE_DEBUG;

    parser.parse_string("output_level", stmp, "info", flag);
    if (flag == 0) driver_params.output_level = get_verbose(stmp);
    librpa::set_output_level(driver_params.output_level);

    _parse_double(opts, vq_threshold);

    _parse_switch(opts, use_kpara_scf_eigvec);

    _parse_string_post(opts, tfgrids_type, get_tfgrid_type);
    if (flag != 0) // backward compatible
    {
        parser.parse_string("tfgrid_type", stmp, flag);
        if (flag == 0) opts.tfgrids_type = get_tfgrid_type(stmp);
    }
    if (opts.tfgrids_type == LIBRPA_TFGRID_UNSET)
        opts.tfgrids_type = LIBRPA_TFGRID_MINIMAX;
    _parse_int(opts, nfreq);
    _parse_double(opts, tfgrids_freq_min);
    _parse_double(opts, tfgrids_freq_interval);
    _parse_double(opts, tfgrids_freq_max);
    _parse_double(opts, tfgrids_time_min);
    _parse_double(opts, tfgrids_time_interval);

    _parse_double(opts, minimax_emin);
    _parse_double(opts, minimax_emax);
    _parse_double(opts, minimax_regulation);

    _parse_switch(opts, use_fullcoul_eps);
    _parse_switch(opts, use_fullcoul_exx);
    _parse_switch(opts, use_fullcoul_wc);
    _parse_switch(opts, use_symmetry_exx);
    _parse_switch(opts, use_symmetry_gw);
    parser.parse_bool("use_symmetry_rpa", btmp, flag);
    opts.use_symmetry_rpa = flag == 0 ? get_switch(btmp) : opts.use_symmetry_gw;
    _parse_switch(opts, output_abacus_gw_gf);

    _parse_int(opts, n_bands_chi0);
    _parse_int(opts, n_bands_sigc);
    _parse_int(opts, option_bvk_remap);

    // chi0 related
    _parse_switch(opts, use_shrink_abfs);
    _parse_switch(opts, use_shrink_chi);

    // RPA specific
    parser.parse_double("gf_R_threshold", opts.gf_threshold, flag); // backward compatible
    _parse_double(opts, gf_threshold);
    _parse_int(opts, libri_chi0_collect_s0_chunk);
    if (opts.libri_chi0_collect_s0_chunk < 0)
        throw std::runtime_error("libri_chi0_collect_s0_chunk must be non-negative");
    {
        double bytes_tmp = 0.0;
        parser.parse_double("libri_chi0_collect_max_bytes", bytes_tmp, flag);
        if (flag == 0)
        {
            if (bytes_tmp < 0.0 ||
                bytes_tmp > static_cast<double>(std::numeric_limits<long long>::max()))
            {
                throw std::runtime_error("libri_chi0_collect_max_bytes must be a non-negative integer byte count");
            }
            opts.libri_chi0_collect_max_bytes = static_cast<long long>(bytes_tmp);
        }
    }
    if (opts.libri_chi0_collect_max_bytes < 0)
        throw std::runtime_error("libri_chi0_collect_max_bytes must be non-negative");
    _parse_double(opts, libri_chi0_threshold_C);
    _parse_double(opts, libri_chi0_threshold_G);
    _parse_switch(opts, use_scalapack_ecrpa);

    // EXX specific
    _parse_double(opts, libri_exx_threshold_C);
    _parse_double(opts, libri_exx_threshold_D);
    _parse_double(opts, libri_exx_threshold_V);

    // GW specific
    _parse_int(opts, n_params_anacon);
    _parse_int(opts, n_params_anacon_resample);
    if (opts.n_params_anacon == 0)
        throw std::runtime_error("n_params_anacon must not be zero");
    if (opts.n_params_anacon_resample == 0)
        throw std::runtime_error("n_params_anacon_resample must not be zero");
    _parse_string_post(opts, anacon_tfgrids_type, get_tfgrid_type);
    _parse_int(opts, anacon_nfreq);
    if (opts.anacon_nfreq == 0)
        throw std::runtime_error("anacon_nfreq must not be zero");
    _parse_double(opts, sqrt_coulomb_threshold);
    _parse_switch(opts, use_scalapack_gw_wc);
    _parse_switch(opts, read_sigc_mat_rf);
    if (flag != 0)  // backward-compatible
    {
        parser.parse_bool("read_sigc", btmp, flag);
        if (flag == 0) opts.read_sigc_mat_rf = get_switch(btmp);
    }
    _parse_switch(opts, use_cholesky_gw_wc);
    _parse_switch(opts, use_gpu_replace_scalapack);
    _parse_switch(opts, use_elpa_sqrt_coulomb);
    _parse_switch(opts, use_hedin_shift);
    _parse_int(opts, istate_ref_hedin_shift);
    _parse_double(opts, libri_g0w0_threshold_C);
    _parse_double(opts, libri_g0w0_threshold_G);
    _parse_double(opts, libri_g0w0_threshold_Wc);
    _parse_switch(opts, replace_w_head);
    _parse_int(opts, option_dielect_func);
    _parse_switch(opts, use_2d_dielectric);
    _parse_switch(opts, output_gw_sigc_ks_kf);
    if (flag != 0)  // backward-compatible
    {
        parser.parse_bool("output_gw_sigc_ks_if", btmp, flag);
        if (flag == 0) opts.output_gw_sigc_ks_kf = get_switch(btmp);
    }
    _parse_int(opts, rpa_headwing_body_start);
    if (opts.rpa_headwing_body_start < 0)
        throw std::runtime_error("rpa_headwing_body_start must be non-negative");
    parser.parse_string("rpa_headwing_mode", stmp, "qavg", flag);
    if (flag == 0 || flag == 1)
    {
        if (stmp != "qavg" && stmp != "head_only")
            throw std::runtime_error("rpa_headwing_mode must be qavg or head_only");
        std::strncpy(opts.rpa_headwing_mode, stmp.c_str(), LIBRPA_MAX_STRLEN);
        opts.rpa_headwing_mode[LIBRPA_MAX_STRLEN - 1] = '\0';
    }
    _parse_switch(opts, output_gw_sigc_ks_mat_kf);
    if (flag != 0)  // backward-compatible
    {
        parser.parse_bool("output_gw_sigc_mat", btmp, flag);
        if (flag == 0) opts.output_gw_sigc_ks_mat_kf = get_switch(btmp);
    }
    _parse_switch(opts, output_exx_ks_mat_k);
    _parse_int(opts, istate_output_mat_start);
    _parse_int(opts, istate_output_mat_end);
    if (opts.istate_output_mat_start < 0)
        throw std::runtime_error("istate_output_mat_start must be non-negative");
    if (opts.istate_output_mat_end >= 0 &&
        opts.istate_output_mat_end <= opts.istate_output_mat_start)
        throw std::runtime_error(
            "istate_output_mat_end must be negative or greater than "
            "istate_output_mat_start");
    _parse_switch(opts, output_gw_sigc_mat_kf);
    _parse_switch(opts, output_gw_sigc_mat_rt);
    _parse_switch(opts, output_gw_sigc_mat_rf);
    _parse_switch(opts, output_wc_rf);
    _parse_switch(opts, output_wc_rf_atom_pair);
    _parse_int(opts, ifreq_output_wc_start);
    _parse_int(opts, ifreq_output_wc_end);
    {  // backward-compatible
        int option_output_Wc_Rf_mat = 0;
        parser.parse_int("option_output_Wc_Rf_mat", option_output_Wc_Rf_mat, flag);
        if (flag == 0)
        {
            if (option_output_Wc_Rf_mat < 0 || option_output_Wc_Rf_mat > 2)
                throw std::runtime_error("option_output_Wc_Rf_mat must be 0, 1, or 2");
            opts.output_wc_rf = get_switch(option_output_Wc_Rf_mat > 0);
            opts.ifreq_output_wc_start = 0;
            opts.ifreq_output_wc_end = option_output_Wc_Rf_mat == 1 ? 1 : -1;
        }
    }
    if (opts.ifreq_output_wc_start < 0)
        throw std::runtime_error("ifreq_output_wc_start must be non-negative");
    if (opts.ifreq_output_wc_end >= 0 &&
        opts.ifreq_output_wc_end <= opts.ifreq_output_wc_start)
        throw std::runtime_error("ifreq_output_wc_end must be negative or greater than ifreq_output_wc_start");

    // QPE solver
    _parse_int(opts, option_qpe_solver);
    _parse_int(opts, qpe_solver_n_iter_max);
    _parse_double(opts, qpe_solver_thres);
    _parse_double(opts, qpe_solver_damp_factor);
    _parse_switch(opts, use_qpe_adaptive_damp);
    _parse_switch(opts, use_qpe_legacy_update);
    _parse_switch(opts, override_qpe_solver_nan);

    // Spectral function
    _parse_double(opts, sf_gf_omega_shift);
    _parse_double(opts, sf_sigc_omega_shift);

    validate_input_parameters();
}

#undef _parse_int
#undef _parse_double
#undef _parse_switch
#undef _parse_string_post
