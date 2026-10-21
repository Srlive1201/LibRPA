#include <algorithm>
#include <memory>
#include <sstream>
#include <vector>

#include "../core/g0w0_task.h"
#include "../core/gw_output.h"
#include "crpa_helper.h"
#include "instance_manager.h"
#include "librpa_compute.h"
#include "librpa_gw_output.h"

namespace
{
struct GwOutputData
{
    const librpa_int::MeanField &mf;
    std::vector<librpa_int::Vector3_Order<double>> kfrac;
    std::vector<int> mapping;
    std::vector<librpa_int::matrix> vxc;
    std::vector<double> vexx;
    std::vector<librpa_int::cplxdb> sigc;

    GwOutputData(const librpa_int::Dataset &dataset, const LibrpaGwOutputInput &input)
        : mf(input.band_path ? dataset.mf_band : dataset.mf)
    {
        if (input.n_spins != mf.get_n_spins() || input.n_kpoints != mf.get_n_kpoints() ||
            input.n_states != mf.get_n_states() || input.n_spins <= 0 || input.n_kpoints <= 0 ||
            input.n_states <= 0 || input.i_state_low < 0 ||
            input.i_state_high <= input.i_state_low || input.i_state_high > input.n_states ||
            input.n_output_kpoints <= 0 || !input.fractional_k || !input.vxc_ha || !input.vexx_ha ||
            !input.sigc_ri ||
            (!input.output_to_input_kpoint && input.n_output_kpoints != input.n_kpoints))
            throw LIBRPA_INVALID_ARGUMENT("Invalid GW output layout or energy buffers");
        for (int output = 0; output < input.n_output_kpoints; ++output)
        {
            const int k =
                input.output_to_input_kpoint ? input.output_to_input_kpoint[output] : output;
            if (k < 0 || k >= input.n_kpoints)
                throw LIBRPA_INVALID_ARGUMENT("GW output k-point index is out of range");
            if (input.output_to_input_kpoint) mapping.push_back(k);
            kfrac.emplace_back(input.fractional_k[3 * output], input.fractional_k[3 * output + 1],
                               input.fractional_k[3 * output + 2]);
        }
        vxc.resize(input.n_spins);
        for (int spin = 0; spin < input.n_spins; ++spin)
        {
            vxc[spin].create(input.n_kpoints, input.n_states);
            for (int k = 0; k < input.n_kpoints; ++k)
                for (int state = 0; state < input.n_states; ++state)
                    vxc[spin](k, state) =
                        input.vxc_ha[(static_cast<std::size_t>(spin) * input.n_kpoints + k) *
                                         input.n_states +
                                     state];
        }
        const auto count = static_cast<std::size_t>(input.n_spins) * input.n_kpoints *
                           (input.i_state_high - input.i_state_low);
        vexx.assign(input.vexx_ha, input.vexx_ha + count);
        sigc.reserve(count);
        for (std::size_t i = 0; i < count; ++i)
            sigc.emplace_back(input.sigc_ri[2 * i], input.sigc_ri[2 * i + 1]);
    }
};

const LibrpaGwOutputInput &require_gw_output(const LibrpaGwOutputInput *input)
{
    if (!input) throw LIBRPA_INVALID_ARGUMENT("Missing GW output input");
    return *input;
}

void emit_gw_text(const std::ostringstream &output, LibrpaGwTextSink sink, void *userdata)
{
    if (!sink) throw LIBRPA_INVALID_ARGUMENT("Missing GW output text sink");
    const auto text = output.str();
    if (sink(userdata, text.data(), text.size()))
        throw LIBRPA_RUNTIME_ERROR("GW output text sink failed");
}
}  // namespace

void librpa_print_gw_band_analysis(LibrpaHandler *h, const LibrpaGwOutputInput *input)
{
    const auto &in = require_gw_output(input);
    const GwOutputData data(*librpa_int::api::get_dataset_instance(h), in);
    librpa_int::print_band_analysis(data.mf, data.kfrac, data.mapping, data.vxc, data.vexx,
                                    data.sigc, in.i_state_low, in.i_state_high - in.i_state_low);
}

void librpa_write_gw_energy_table(LibrpaHandler *h, const LibrpaGwOutputInput *input,
                                  LibrpaGwTextSink sink, void *userdata)
{
    const auto &in = require_gw_output(input);
    if (in.n_kpoint_weights != in.n_kpoints || !in.kpoint_weights)
        throw LIBRPA_INVALID_ARGUMENT("GW output k-point weight count does not match loaded data");
    const GwOutputData data(*librpa_int::api::get_dataset_instance(h), in);
    const std::vector<double> weights(in.kpoint_weights, in.kpoint_weights + in.n_kpoint_weights);
    std::ostringstream output;
    librpa_int::write_energy_qp(output, data.mf, data.kfrac, data.mapping, data.vxc, data.vexx,
                                data.sigc, in.n_kpoints, in.i_state_low,
                                in.i_state_high - in.i_state_low, weights);
    emit_gw_text(output, sink, userdata);
}

void librpa_write_gw_band_spin_table(LibrpaHandler *h, const LibrpaGwOutputInput *input,
                                     const int spin, LibrpaGwTextSink sink, void *userdata)
{
    const auto &in = require_gw_output(input);
    if (in.output_to_input_kpoint)
        throw LIBRPA_INVALID_ARGUMENT("GW band spin table requires the loaded k-point order");
    const GwOutputData data(*librpa_int::api::get_dataset_instance(h), in);
    std::ostringstream output;
    librpa_int::write_gw_band_spin_table(output, data.mf, data.kfrac, data.vxc, data.vexx,
                                         data.sigc, spin, in.i_state_low,
                                         in.i_state_high - in.i_state_low);
    emit_gw_text(output, sink, userdata);
}

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
    librpa_int::build_g0w0_sigma(
        *dataset, *options,
        [&](const librpa_int::Chi0QMap &response)
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
    librpa_int::get_g0w0_qpe_kgrid(*dataset, *p_opts, n_spins, n_kpts_this, iks_this, i_state_low,
                                   i_state_high, vxc, vexx, sigc_re, sigc_im, eqp);
}

void librpa_get_g0w0_spectral_function_kgrid(
    LibrpaHandler *h, const LibrpaOptions *p_opts, const int n_spins,
    const int n_kpts_this, const int *iks_this, int i_state_low, int i_state_high,
    const int n_omegas, const double *omegas, const double *vxc, const double *vexx,
    double *spectral_function, double *sigc)
{
    if (!p_opts) throw LIBRPA_INVALID_ARGUMENT("missing GW options");
    const auto dataset = librpa_int::api::get_dataset_instance(h);
    librpa_int::get_g0w0_spectral_function_kgrid(*dataset, *p_opts, n_spins, n_kpts_this, iks_this,
                                                 i_state_low, i_state_high, n_omegas, omegas, vxc,
                                                 vexx, spectral_function, sigc);
}

void librpa_get_g0w0_qpe_band_k(LibrpaHandler *h, const LibrpaOptions *p_opts, const int n_spins,
                                const int n_kpts_band_this, const int *iks_band_this,
                                int i_state_low, int i_state_high, const double *vxc_band,
                                const double *vexx_band, double *sigc_band_re,
                                double *sigc_band_im, double *eqp_band)
{
    if (!p_opts) throw LIBRPA_INVALID_ARGUMENT("missing GW options");
    const auto dataset = librpa_int::api::get_dataset_instance(h);
    librpa_int::get_g0w0_qpe_band_k(*dataset, *p_opts, n_spins, n_kpts_band_this, iks_band_this,
                                    i_state_low, i_state_high, vxc_band, vexx_band, sigc_band_re,
                                    sigc_band_im, eqp_band);
}

void librpa_get_g0w0_spectral_function_band_k(LibrpaHandler *h, const LibrpaOptions *p_opts,
                                              const int n_spins, const int n_kpts_band_this,
                                              const int *iks_band_this, int i_state_low,
                                              int i_state_high, const int n_omegas,
                                              const double *omegas, const double *vxc_band,
                                              const double *vexx_band,
                                              double *spectral_function_band, double *sigc_band)
{
    if (!p_opts) throw LIBRPA_INVALID_ARGUMENT("missing GW options");
    const auto dataset = librpa_int::api::get_dataset_instance(h);
    librpa_int::get_g0w0_spectral_function_band_k(
        *dataset, *p_opts, n_spins, n_kpts_band_this, iks_band_this, i_state_low, i_state_high,
        n_omegas, omegas, vxc_band, vexx_band, spectral_function_band, sigc_band);
}

void librpa_get_g0w0_sigc_kgrid(LibrpaHandler *h, const LibrpaOptions *p_opts, const int n_spins,
                                const int n_kpts_this, const int *iks_this, int i_state_low,
                                int i_state_high, const double *vxc, const double *vexx,
                                double *sigc_re, double *sigc_im)
{
    const int n_states_calc = std::max(0, i_state_high - i_state_low);
    std::vector<double> eqp(n_spins * n_kpts_this * n_states_calc);
    librpa_get_g0w0_qpe_kgrid(h, p_opts, n_spins, n_kpts_this, iks_this, i_state_low, i_state_high,
                              vxc, vexx, sigc_re, sigc_im, eqp.data());
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
