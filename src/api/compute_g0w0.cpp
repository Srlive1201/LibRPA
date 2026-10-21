#include "librpa_compute.h"

#include <algorithm>
#include <memory>
#include <vector>

#include "../core/g0w0_task.h"
#include "crpa_helper.h"
#include "instance_manager.h"

void librpa_build_g0w0_sigma(LibrpaHandler *h, const LibrpaOptions *options)
{
    if (!options) throw LIBRPA_INVALID_ARGUMENT("missing GW options");
    const auto dataset = librpa_int::api::get_dataset_instance(h);
    librpa_int::build_g0w0_sigma(*dataset, *options);
}

LibrpaCrpaResult *librpa_build_g0w0_sigma_crpa(LibrpaHandler *h, const LibrpaOptions *options,
                                               const LibrpaCrpaInput *input)
{
    if (!options || !input) throw LIBRPA_INVALID_ARGUMENT("missing shared GW/cRPA input");
    const auto dataset = librpa_int::api::get_dataset_instance(h);
    auto model = librpa_int::prepare_crpa_task_model(
        *dataset, *options, librpa_int::api::unpack_crpa_input(*dataset, input),
        librpa_int::CrpaTaskMode::SharedGw);
    std::unique_ptr<LibrpaCrpaResult, decltype(&librpa_delete_crpa_result)> result(
        nullptr, librpa_delete_crpa_result);
    librpa_int::build_g0w0_sigma(*dataset, *options, [&](const librpa_int::Chi0QMap &response)
    {
        result.reset(librpa_int::api::own_crpa_result(
            librpa_int::evaluate_crpa_task_from_response(*dataset, *options, model, response)));
    });
    return result.release();
}

void librpa_get_g0w0_qpe_kgrid(LibrpaHandler *h, const LibrpaOptions *p_opts, const int n_spins,
                               const int n_kpts_this, const int *iks_this, int i_state_low,
                               int i_state_high, const double *vxc, const double *vexx,
                               double *sigc_re, double *sigc_im, double *eqp)
{
    if (!p_opts) throw LIBRPA_INVALID_ARGUMENT("missing GW options");
    const auto dataset = librpa_int::api::get_dataset_instance(h);
    librpa_int::get_g0w0_qpe_kgrid(*dataset, *p_opts, n_spins, n_kpts_this, iks_this, i_state_low, i_state_high, vxc, vexx, sigc_re, sigc_im, eqp);
}

void librpa_get_g0w0_spectral_function_kgrid(
    LibrpaHandler *h, const LibrpaOptions *p_opts, const int n_spins,
    const int n_kpts_this, const int *iks_this, int i_state_low, int i_state_high,
    const int n_omegas, const double *omegas, const double *vxc, const double *vexx,
    double *spectral_function, double *sigc)
{
    if (!p_opts) throw LIBRPA_INVALID_ARGUMENT("missing GW options");
    const auto dataset = librpa_int::api::get_dataset_instance(h);
    librpa_int::get_g0w0_spectral_function_kgrid(*dataset, *p_opts, n_spins, n_kpts_this, iks_this, i_state_low, i_state_high, n_omegas, omegas, vxc, vexx, spectral_function, sigc);
}

void librpa_get_g0w0_qpe_band_k(LibrpaHandler *h, const LibrpaOptions *p_opts, const int n_spins,
                                const int n_kpts_band_this, const int *iks_band_this,
                                int i_state_low, int i_state_high, const double *vxc_band,
                                const double *vexx_band, double *sigc_band_re,
                                double *sigc_band_im, double *eqp_band)
{
    if (!p_opts) throw LIBRPA_INVALID_ARGUMENT("missing GW options");
    const auto dataset = librpa_int::api::get_dataset_instance(h);
    librpa_int::get_g0w0_qpe_band_k(*dataset, *p_opts, n_spins, n_kpts_band_this, iks_band_this, i_state_low, i_state_high, vxc_band, vexx_band, sigc_band_re, sigc_band_im, eqp_band);
}

void librpa_get_g0w0_spectral_function_band_k(
    LibrpaHandler *h, const LibrpaOptions *p_opts, const int n_spins,
    const int n_kpts_band_this, const int *iks_band_this, int i_state_low,
    int i_state_high, const int n_omegas, const double *omegas,
    const double *vxc_band, const double *vexx_band, double *spectral_function_band,
    double *sigc_band)
{
    if (!p_opts) throw LIBRPA_INVALID_ARGUMENT("missing GW options");
    const auto dataset = librpa_int::api::get_dataset_instance(h);
    librpa_int::get_g0w0_spectral_function_band_k(*dataset, *p_opts, n_spins, n_kpts_band_this, iks_band_this, i_state_low, i_state_high, n_omegas, omegas, vxc_band, vexx_band, spectral_function_band, sigc_band);
}

void librpa_get_g0w0_sigc_kgrid(LibrpaHandler *h, const LibrpaOptions *p_opts, const int n_spins,
                                const int n_kpts_this, const int *iks_this, int i_state_low,
                                int i_state_high, const double *vxc, const double *vexx,
                                double *sigc_re, double *sigc_im)
{
    const int n_states_calc = std::max(0, i_state_high - i_state_low);
    std::vector<double> eqp(n_spins * n_kpts_this * n_states_calc);
    librpa_get_g0w0_qpe_kgrid(h, p_opts, n_spins, n_kpts_this, iks_this, i_state_low,
                              i_state_high, vxc, vexx, sigc_re, sigc_im, eqp.data());
}

void librpa_get_g0w0_sigc_band_k(LibrpaHandler *h, const LibrpaOptions *p_opts, const int n_spins,
                                 const int n_kpts_band_this, const int *iks_band_this,
                                 int i_state_low, int i_state_high, const double *vxc_band,
                                 const double *vexx_band, double *sigc_band_re,
                                 double *sigc_band_im)
{
    const int n_states_calc = std::max(0, i_state_high - i_state_low);
    std::vector<double> eqp_band(n_spins * n_kpts_band_this * n_states_calc);
    librpa_get_g0w0_qpe_band_k(h, p_opts, n_spins, n_kpts_band_this, iks_band_this,
                               i_state_low, i_state_high, vxc_band, vexx_band, sigc_band_re,
                               sigc_band_im, eqp_band.data());
}
