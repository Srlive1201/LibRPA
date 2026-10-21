#include <librpa_enums.h>

#include <algorithm>
#include <cmath>
#include <exception>
#include <fstream>
#include <iomanip>
#include <map>
#include <numeric>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "../../src/api/compute_helper.h"
#include "../../src/api/dataset_helper.h"
#include "../../src/api/instance_manager.h"
#include "../../src/core/qsgw/convergence.h"
#include "../../src/core/qsgw/correlation_potential.h"
#include "../../src/core/qsgw/distributed_matrix.h"
#include "../../src/core/qsgw/effective_hamiltonian.h"
#include "../../src/core/qsgw/fixed_basis.h"
#include "../../src/core/qsgw/hamiltonian_cut.h"
#include "../../src/core/qsgw/hamiltonian_mixing.h"
#include "../../src/core/qsgw/headwing_velocity.h"
#include "../../src/core/qsgw/occupation.h"
#include "../../src/core/qsgw/projection_target.h"
#include "../../src/io/fs.h"
#include "../../src/io/global_io.h"
#include "../../src/math/utils_matrix_mpi.h"
#include "../../src/mpi/collective_stage.h"
#include "../../src/utils/constants.h"
#include "../../src/utils/profiler.h"
#include "../driver.h"
#include "../qsgw/ao_operator_io.h"
#include "../qsgw/band_output.h"
#include "../qsgw/checkpoint.h"
#include "../qsgw/iteration_trace.h"
#include "../qsgw/vxc_io.h"
#include "../read_data.h"
#include "../task.h"
#include "task_helper.h"

namespace
{

using librpa_int::Matz;
using librpa_int::MeanField;
using librpa_int::Vector3_Order;
using librpa_int::cplxdb;
using librpa_int::qsgw::ScopedReferenceEigenvectors;
using librpa_int::qsgw::SpinKMatrixMap;
using librpa_int::qsgw::VelocityMatrix;
using SigmaMatrixMap = librpa_int::qsgw::SpinKFrequencyMatrixMap;

class ScopedSigmaMatrixRetention
{
public:
    explicit ScopedSigmaMatrixRetention(bool& flag)
        : flag_(flag), original_(flag_)
    {
        flag_ = true;
    }

    ~ScopedSigmaMatrixRetention() noexcept
    {
        flag_ = original_;
    }

    ScopedSigmaMatrixRetention(const ScopedSigmaMatrixRetention&) = delete;
    ScopedSigmaMatrixRetention& operator=(
        const ScopedSigmaMatrixRetention&) = delete;

private:
    bool& flag_;
    bool original_;
};

class ScopedLegacyHeadOccupations
{
public:
    explicit ScopedLegacyHeadOccupations(MeanField& meanfield)
        : meanfield_(meanfield),
          live_weights_(meanfield.get_weight())
    {
        auto& weights = meanfield_.get_weight();
        for (const auto& spin_weights : weights)
        {
            if (spin_weights.nr != meanfield_.get_n_kpoints() ||
                spin_weights.nc != meanfield_.get_n_bands())
            {
                throw std::invalid_argument(
                    "QSGW head occupation matrix has an invalid shape");
            }
        }
        for (auto& spin_weights : weights)
        {
            for (int kpoint = 1; kpoint < spin_weights.nr; ++kpoint)
            {
                for (int band = 0; band < spin_weights.nc; ++band)
                    spin_weights(kpoint, band) = spin_weights(0, band);
            }
        }
    }

    ~ScopedLegacyHeadOccupations() noexcept
    {
        meanfield_.get_weight().swap(live_weights_);
    }

    ScopedLegacyHeadOccupations(
        const ScopedLegacyHeadOccupations&) = delete;
    ScopedLegacyHeadOccupations& operator=(
        const ScopedLegacyHeadOccupations&) = delete;

private:
    MeanField& meanfield_;
    std::vector<librpa_int::matrix> live_weights_;
};

template <typename Function>
void collective_root_stage(
    const librpa_int::MpiCommHandler& communicator,
    const std::string& label,
    Function&& function)
{
    std::exception_ptr root_error;
    int failed = 0;
    if (communicator.is_root())
    {
        try
        {
            function();
        }
        catch (const std::exception& error)
        {
            failed = 1;
            root_error = std::current_exception();
            librpa_int::global::lib_printf(
                LIBRPA_VERBOSE_CRITICAL, "%s failed: %s\n",
                label.c_str(), error.what());
        }
        catch (...)
        {
            failed = 1;
            root_error = std::current_exception();
        }
    }
    communicator.bcast(&failed, 1, 0);
    if (failed)
    {
        if (communicator.is_root()) std::rethrow_exception(root_error);
        throw LIBRPA_RUNTIME_ERROR(label + " failed on root");
    }
}

std::string resolve_input_path(const std::string& base,
                               const std::string& path)
{
    return librpa_int::is_absolute_path(path)
               ? path
               : librpa_int::join_path(base, path);
}

void prepare_qsgw_symmetry_context(librpa_int::Dataset& dataset)
{
    const bool symmetry_reduced_scf_grid =
        dataset.pbc.kfrac_list.size() < dataset.pbc.kfrac_list_full.size();
    if (!symmetry_reduced_scf_grid) return;

    const bool all_symmetry_routes_enabled =
        driver::get_bool(driver::opts.use_symmetry_exx) &&
        driver::get_bool(driver::opts.use_symmetry_gw) &&
        driver::get_bool(driver::opts.use_symmetry_rpa);
    if (all_symmetry_routes_enabled)
        librpa_int::initialize_symmetry_context(dataset, true);
}

void validate_qsgw_reference(const librpa_int::Dataset& dataset,
                             const bool compute_band)
{
    using namespace librpa_int::qsgw;
    validate_projection_target(
        dataset.mf, dataset.pbc.kfrac_list,
        dataset.mf.get_n_spins(), dataset.mf.get_n_spinor(),
        dataset.mf.get_n_aos(), "grid");
    const bool symmetry_reduced_scf_grid =
        dataset.pbc.kfrac_list.size() < dataset.pbc.kfrac_list_full.size();
    if (symmetry_reduced_scf_grid)
    {
        const bool all_symmetry_routes_enabled =
            driver::get_bool(driver::opts.use_symmetry_exx) &&
            driver::get_bool(driver::opts.use_symmetry_gw) &&
            driver::get_bool(driver::opts.use_symmetry_rpa);
        if (!all_symmetry_routes_enabled ||
            !dataset.symmetry_context.available ||
            dataset.symmetry_context.kstars.size() !=
                dataset.pbc.kfrac_list.size() ||
            dataset.symmetry_context.count_kstar_members() !=
                dataset.pbc.kfrac_list_full.size())
        {
            throw std::invalid_argument(
                "QSGW symmetry-reduced SCF input requires complete EXX/GW/RPA k-star restoration");
        }
    }
    if (compute_band)
    {
        const int complete_dimension =
            dataset.mf.get_n_aos() * dataset.mf.get_n_spinor();
        if (dataset.mf.get_n_bands() != complete_dimension ||
            dataset.mf_band.get_n_bands() != complete_dimension)
        {
            throw std::invalid_argument(
                "QSGW fixed-basis band rotation requires complete square grid and band references");
        }
        validate_projection_target(
            dataset.mf_band, dataset.kfrac_band_list,
            dataset.mf.get_n_spins(), dataset.mf.get_n_spinor(),
            dataset.mf.get_n_aos(), "band");
    }
}

SpinKMatrixMap copy_exx_root(const librpa_int::Exx& exchange,
                             const MeanField& reference)
{
    SpinKMatrixMap result;
    for (int spin = 0; spin < reference.get_n_spins(); ++spin)
    {
        for (int kpoint = 0; kpoint < reference.get_n_kpoints(); ++kpoint)
        {
            const auto spin_it = exchange.exx_KS.find(spin);
            if (spin_it == exchange.exx_KS.end())
                throw std::runtime_error("QSGW EXX spin channel is missing");
            const auto k_it = spin_it->second.find(kpoint);
            if (k_it == spin_it->second.end() ||
                k_it->second.nr() != reference.get_n_bands() ||
                k_it->second.nc() != reference.get_n_bands())
            {
                throw std::runtime_error(
                    "QSGW EXX fixed-basis matrix is missing or malformed");
            }
            result[spin][kpoint] = k_it->second.copy();
        }
    }
    return result;
}

SigmaMatrixMap collect_sigma_root(
    const librpa_int::G0W0& self_energy,
    const MeanField& reference,
    const std::vector<double>& frequencies,
    const librpa_int::MpiCommHandler& communicator)
{
    using librpa_int::collect_blacs_matrix_root;
    SigmaMatrixMap result;
    for (int spin = 0; spin < reference.get_n_spins(); ++spin)
    {
        for (int kpoint = 0; kpoint < reference.get_n_kpoints(); ++kpoint)
        {
            for (const double frequency : frequencies)
            {
                const Matz* local = nullptr;
                const auto spin_it = self_energy.sigc_is_ik_f_KS.find(spin);
                if (spin_it != self_energy.sigc_is_ik_f_KS.end())
                {
                    const auto k_it = spin_it->second.find(kpoint);
                    if (k_it != spin_it->second.end())
                    {
                        const auto f_it = k_it->second.find(frequency);
                        if (f_it != k_it->second.end()) local = &f_it->second;
                    }
                }
                int valid = local != nullptr &&
                                    local->major() == librpa_int::MAJOR::COL &&
                                    local->nr() ==
                                        self_energy.desc_sigc_is_ik_f_KS.m_loc() &&
                                    local->nc() ==
                                        self_energy.desc_sigc_is_ik_f_KS.n_loc()
                                ? 1
                                : 0;
                int all_valid = 0;
                communicator.allreduce(&valid, &all_valid, 1, MPI_MIN);
                if (!all_valid)
                    throw std::runtime_error(
                        "QSGW distributed fixed-basis Sigma is incomplete");
                Matz full = collect_blacs_matrix_root(
                    *local, self_energy.desc_sigc_is_ik_f_KS);
                if (communicator.is_root())
                {
                    if (full.nr() != reference.get_n_bands() ||
                        full.nc() != reference.get_n_bands())
                    {
                        throw std::runtime_error(
                            "QSGW collected Sigma has an invalid shape");
                    }
                    result[spin][kpoint][frequency] = std::move(full);
                }
            }
        }
    }
    return result;
}

SpinKMatrixMap build_correlation_map(
    const MeanField& live,
    const SigmaMatrixMap& sigma,
    const std::vector<double>& frequencies,
    const LibrpaOptions& options)
{
    using namespace librpa_int::qsgw;
    CorrelationPotentialSettings settings;
    settings.mode = CorrelationPotentialMode::ModeB;
    settings.n_params_anacon = options.n_params_anacon;
    settings.n_params_anacon_resample =
        options.n_params_anacon_resample < 0
            ? options.n_params_anacon
            : options.n_params_anacon_resample;

    SpinKMatrixMap result;
    for (int spin = 0; spin < live.get_n_spins(); ++spin)
    {
        for (int kpoint = 0; kpoint < live.get_n_kpoints(); ++kpoint)
        {
            result[spin][kpoint] = build_qsgw_correlation_potential(
                live, frequencies, sigma.at(spin).at(kpoint),
                spin, kpoint, settings);
        }
    }
    return result;
}

// The matrices are already available for QSGW. Rotate to the entering
// eigenframe and reuse G0W0 continuation/roots without a second GW build.
void write_qsgw_qp_state(std::ostream& output, int state_index, int update,
                         const std::string& channel, const MeanField& live,
                         const std::vector<Vector3_Order<double>>& kpoints,
                         const SpinKMatrixMap& rotation, const SpinKMatrixMap& entering_hamiltonian,
                         const SpinKMatrixMap& reference_hamiltonian, const SpinKMatrixMap& dft_vxc,
                         const SpinKMatrixMap& exchange, const SigmaMatrixMap& sigma,
                         const librpa_int::TFGrids& tfg, const std::vector<double>& kpoint_weights)
{
    const auto& options = driver::opts;
    const int low = std::max(0, driver::driver_params.i_state_low);
    const int high = std::min(live.get_n_states(), driver::driver_params.i_state_high);
    const int n = high - low, nk = live.get_n_kpoints(), ns = live.get_n_spins();
    if (n <= 0) return;
    const auto frequencies = tfg.get_freq_nodes();
    std::vector<int> iks(nk);
    std::iota(iks.begin(), iks.end(), 0);
    std::vector<double> potential(ns * nk * n), exx(potential.size());
    std::vector<cplxdb> diagonal;
    diagonal.reserve(potential.size() * frequencies.size());
    for (int spin = 0; spin < ns; ++spin)
        for (int kpoint = 0; kpoint < nk; ++kpoint)
        {
            const auto& u = rotation.at(spin).at(kpoint);
            const auto adjoint = transpose(u, true);
            const auto rotate = [&](const Matz& matrix)
            {
                // BLACS Sigma/EXX are column-major; the fixed-frame rotation
                // can be row-major. Normalize storage before the BLAS product.
                const Matz aligned(matrix.nr(), matrix.nc(), matrix.ptr(), matrix.major(),
                                   u.major());
                return adjoint * aligned * u;
            };
            // H_t = H_DFT - Vxc_DFT + Vstatic_t, including mixing/cuts/restart.
            Matz static_reference(u.nr(), u.nc(), u.major());
            for (int row = 0; row < u.nr(); ++row)
                for (int column = 0; column < u.nc(); ++column)
                    static_reference(row, column) =
                        entering_hamiltonian.at(spin).at(kpoint)(row, column) -
                        reference_hamiltonian.at(spin).at(kpoint)(row, column) +
                        dft_vxc.at(spin).at(kpoint)(row, column);
            const auto vstatic = rotate(static_reference);
            const auto sigx = rotate(exchange.at(spin).at(kpoint));
            for (int band = low; band < high; ++band)
            {
                const int i = (spin * nk + kpoint) * n + band - low;
                potential[i] = vstatic(band, band).real();
                exx[i] = sigx(band, band).real();
            }
            for (double frequency : frequencies)
            {
                const auto sigc = rotate(sigma.at(spin).at(kpoint).at(frequency));
                for (int band = low; band < high; ++band) diagonal.push_back(sigc(band, band));
            }
        }
    std::vector<double> real(potential.size()), imag(potential.size()), eqp(potential.size());
    std::vector<int> info(potential.size());
    librpa_int::api::solve_g0w0_qpe_diagonal(
        live, tfg, options, diagonal, iks, nk, iks.data(), low, n, potential.data(), exx.data(),
        real.data(), imag.data(), eqp.data(), "QSGW entering-state G0W0", info.data());
    for (int spin = 0; spin < ns; ++spin)
        for (int kpoint = 0; kpoint < nk; ++kpoint)
            for (int band = low; band < high; ++band)
            {
                const int i = (spin * nk + kpoint) * n + band - low;
                const auto& k = kpoints.at(kpoint);
                output << state_index << ',' << update << ',' << channel << ',' << spin + 1 << ','
                       << kpoint + 1 << ',' << k.x << ',' << k.y << ',' << k.z << ',' << band + 1
                       << ',' << live.get_weight().at(spin)(kpoint, band) * nk << ','
                       << live.get_eigenvals().at(spin)(kpoint, band) * librpa_int::HA2EV << ','
                       << potential[i] * librpa_int::HA2EV << ',' << exx[i] * librpa_int::HA2EV
                       << ',' << real[i] * librpa_int::HA2EV << ',' << imag[i] * librpa_int::HA2EV
                       << ',' << eqp[i] * librpa_int::HA2EV << ',' << info[i] << '\n';
            }
    output.flush();
    if (!output) throw LIBRPA_RUNTIME_ERROR("Cannot write QSGW quasiparticle output");
    // Use the ordinary G0W0 writers, including their energy/occupation units
    // and selected-state layout. The iteration is encoded in the filename.
    std::vector<librpa_int::matrix> vstatic_diagonal(ns);
    std::vector<cplxdb> correlation(potential.size());
    for (int spin = 0; spin < ns; ++spin)
    {
        vstatic_diagonal[spin].create(nk, live.get_n_states());
        vstatic_diagonal[spin].zero_out();
        for (int kpoint = 0; kpoint < nk; ++kpoint)
            for (int band = low; band < high; ++band)
            {
                const int i = (spin * nk + kpoint) * n + band - low;
                vstatic_diagonal[spin](kpoint, band) = potential[i];
                correlation[i] = cplxdb(real[i], imag[i]);
            }
    }
    if (channel == "grid")
        write_energy_qp(
            live, kpoints, {}, vstatic_diagonal, exx, correlation, nk, low, n, kpoint_weights,
            librpa_int::join_path(options.output_dir, "energy_qp_" + std::to_string(update)));
    if (channel == "band")
        for (int spin = 0; spin < ns; ++spin)
        {
            std::ostringstream name;
            name << "GW_band_spin_" << spin + 1 << '_' << update << ".dat";
            std::ofstream band_output(librpa_int::join_path(options.output_dir, name.str()));
            if (!band_output)
                throw LIBRPA_RUNTIME_ERROR("Cannot open QSGW iteration GW band output");
            write_gw_band_spin_table(band_output, live, kpoints, vstatic_diagonal, exx, correlation,
                                     spin, low, n);
        }
}

void run_qsgw()
{
    using namespace driver;
    using namespace librpa_int;
    using namespace librpa_int::global;
    using namespace librpa_int::qsgw;

    const bool update_head = driver::get_bool(opts.replace_w_head) &&
                             (opts.option_dielect_func == 3 || opts.option_dielect_func == 4);
    const bool update_wing = update_head && opts.option_dielect_func == 3;
    if (driver_params.use_pyatb)
    {
        throw LIBRPA_RUNTIME_ERROR(
            "QSGW independent PyATB head updates are unsupported; use the same-grid velocity input");
    }
    if (driver::get_bool(opts.replace_w_head) && !update_head)
    {
        throw LIBRPA_RUNTIME_ERROR("QSGW analytic head/wing requires option_dielect_func = 3 or 4");
    }
    if (driver::get_bool(opts.use_kpara_scf_eigvec))
    {
        throw LIBRPA_RUNTIME_ERROR(
            "QSGW fixed-basis iteration requires replicated SCF wavefunctions; set use_kpara_scf_eigvec = false");
    }

    profiler.start("qsgw", "QSGW fixed-basis self-consistent calculation");
    const auto dataset = api::get_dataset_instance(h);
    if (update_head && driver_params.qsgw_headwing_velocity == "current_hr")
    {
        const auto& mf = dataset->mf;
        if ((mf.get_n_spinor() != 1 && mf.get_n_spinor() != 2) ||
            (mf.get_n_spinor() == 2 && mf.get_n_spins() != 1) ||
            mf.get_n_bands() != mf.get_n_aos() * mf.get_n_spinor())
            throw LIBRPA_RUNTIME_ERROR(
                "QSGW current-HR velocity requires a complete scalar or two-spinor AO basis; "
                "SOC has one physical spin channel");
        if (dataset->pbc.kfrac_list.size() != dataset->pbc.kfrac_list_full.size() ||
            mf.get_n_kpoints() != dataset->pbc.get_n_cells_bvk())
            throw LIBRPA_RUNTIME_ERROR(
                "QSGW current-HR velocity requires the complete full-BZ SCF mesh; "
                "IBZ velocity sewing is unsupported");
        if (mf.get_n_spinor() == 2 &&
            (driver::get_bool(opts.use_symmetry_exx) || driver::get_bool(opts.use_symmetry_gw) ||
             driver::get_bool(opts.use_symmetry_rpa)))
            throw LIBRPA_RUNTIME_ERROR(
                "QSGW SOC current-HR velocity requires EXX/GW/RPA symmetry switches off; "
                "magnetic-symmetry velocity sewing is unsupported");
    }
    if (opts.parallel_routing != LIBRPA_ROUTING_LIBRI &&
        opts.parallel_routing != LIBRPA_ROUTING_AUTO)
    {
        throw LIBRPA_RUNTIME_ERROR("QSGW requires parallel_routing=libri or auto");
    }
    const std::string band_kpath_path =
        resolve_input_path(driver_params.input_dir, driver_params.fn_band_kpath_info);
    int has_band_kpath = 0;
    collective_root_stage(
        dataset->comm_h, "QSGW band-path discovery",
        [&]
        {
            has_band_kpath =
                !driver_params.fn_band_kpath_info.empty() && path_exists(band_kpath_path.c_str());
            if (has_band_kpath)
                lib_printf("QSGW band k-path file %s found; include band-path updates\n",
                           band_kpath_path.c_str());
            else
                lib_printf("QSGW band k-path file %s not found; use SCF-grid updates only\n",
                           band_kpath_path.c_str());
        });
    dataset->comm_h.bcast(&has_band_kpath, 1, 0);
    const bool compute_band = has_band_kpath != 0;
    if (compute_band)
    {
        collective_all_stage(dataset->comm_h, "QSGW band k-path input",
                             [&] { read_band_kpath_info(band_kpath_path); });
        dataset->comm_h.barrier();
        collective_all_stage(dataset->comm_h, "QSGW band mean-field input",
                             [&] { read_band_meanfield_data(driver_params.input_dir); });
        dataset->comm_h.barrier();
    }
    read_Vq_row(driver_params.input_dir, driver_params.prefix_coul_cut,
                opts.vq_threshold, local_atpair, true,
                driver_params.version_coul_reader,
                driver::get_bool(opts.use_shrink_abfs));
    collective_root_stage(dataset->comm_h, "QSGW reference validation",
                          [&]
                          {
                              prepare_qsgw_symmetry_context(*dataset);
                              validate_qsgw_reference(*dataset, compute_band);
                          });

    const MeanField reference = dataset->mf;
    // Use the full SCF mesh as well as the loaded count: one IBZ representative
    // of a multi-k mesh must not enter the single-k supercell branch.
    const bool single_kpoint =
        reference.get_n_kpoints() == 1 && dataset->pbc.get_n_cells_bvk() == 1;
    const std::string effective_mixer = single_kpoint ? driver_params.qsgw_mixer : "none";
    const double mixing_beta = effective_mixer == "linear" ? driver_params.qsgw_mixing_beta : 1.0;
    if (dataset->comm_h.is_root())
        lib_printf(
            "QSGW mixing from SCF input: loaded_kpoints=%d mesh=%dx%dx%d requested=%s effective=%s "
            "beta=%.8f\n",
            reference.get_n_kpoints(), dataset->pbc.period.x, dataset->pbc.period.y,
            dataset->pbc.period.z, driver_params.qsgw_mixer.c_str(), effective_mixer.c_str(),
            mixing_beta);
    const double electron_count = physical_electron_count(
        reference, dataset->pbc.weight_k);
    const OccupationResult initial_occupations = analyze_qsgw_occupations(
        reference, dataset->pbc.weight_k, electron_count);
    const SpinKMatrixMap reference_hamiltonian =
        build_reference_hamiltonian(reference);

    VelocityMatrix reference_velocity;
    std::optional<HeadwingAoOperators> ao_operators;
    const bool current_hr_velocity =
        update_head && driver_params.qsgw_headwing_velocity == "current_hr";
    if (current_hr_velocity)
    {
        collective_all_stage(
            dataset->comm_h, "QSGW head/wing AO input",
            [&]
            {
                ao_operators = driver::qsgw::read_headwing_ao_operators(
                    resolve_input_path(driver_params.input_dir, driver_params.qsgw_fn_sr),
                    resolve_input_path(driver_params.input_dir, driver_params.qsgw_fn_rr));
            });
    }
    else if (driver_params.qsgw_export_hamiltonian_for_pyatb)
    {
        collective_root_stage(
            dataset->comm_h, "QSGW HR export AO input",
            [&]
            {
                ao_operators = driver::qsgw::read_headwing_ao_operators(
                    resolve_input_path(driver_params.input_dir, driver_params.qsgw_fn_sr),
                    resolve_input_path(driver_params.input_dir, driver_params.qsgw_fn_rr));
            });
    }
    if (current_hr_velocity || driver_params.qsgw_export_hamiltonian_for_pyatb)
    {
        const auto load_native_hamiltonian = [&]
        {
            // Export and current_hr differentiate the same native-baseline
            // reconstruction of the complete entering/current operator.
            for (int spin = 0; spin < reference.get_n_spins(); ++spin)
            {
                const std::string filename = "hrs" + std::to_string(spin + 1) + "_nao.csr";
                ao_operators->native_dft_hamiltonian.push_back(
                    driver::qsgw::read_headwing_hamiltonian_csr(
                        resolve_input_path(driver_params.input_dir, filename), 0.5));
            }
            if (dataset->comm_h.is_root())
                lib_printf(
                    "QSGW native-baseline H(R): preserve DFT H(R), Fourier-transform "
                    "only the complete current correction for HR export/current_hr; "
                    "%d physical spin channel(s)\n",
                    reference.get_n_spins());
        };
        if (current_hr_velocity)
            collective_all_stage(dataset->comm_h, "QSGW native DFT Hamiltonian input",
                                 load_native_hamiltonian);
        else
            collective_root_stage(dataset->comm_h, "QSGW native DFT Hamiltonian input",
                                  load_native_hamiltonian);
    }
    if (update_head && driver_params.qsgw_headwing_velocity != "current_hr")
    {
        read_headwing_input(driver_params.input_dir, false);
        reference_velocity = dataset->velocity_matrix;
        if (driver_params.constants_choice == "aims")
        {
            // FHI-aims head input uses only interband velocity elements.
            prepare_interband_velocity(reference_velocity, reference);
            align_distributed_velocity_to_reference_wfc(dataset->p_headwing->get_meanfield_df(),
                                                        reference, reference_velocity,
                                                        dataset->comm_h);
        }
    }

    std::optional<MeanField> band_reference;
    SpinKMatrixMap band_reference_hamiltonian;
    if (compute_band)
    {
        dataset->mf_band.get_efermi() = reference.get_efermi();
        band_reference = dataset->mf_band;
        band_reference_hamiltonian =
            build_reference_hamiltonian(*band_reference);
    }

    SpinKMatrixMap dft_vxc;
    SpinKMatrixMap dft_vxc_band;
    collective_root_stage(dataset->comm_h, "QSGW Vxc input", [&] {
        const bool aims_input = driver_params.constants_choice == "aims";
        const VxcBasis basis = driver_params.qsgw_vxc_basis == "nao"
            ? VxcBasis::Nao : VxcBasis::State;
        dft_vxc = read_qsgw_vxc(
            driver_params.input_dir, driver_params.prefix_vxc_scf,
            reference, aims_input, false, basis);
        if (compute_band)
            dft_vxc_band = read_qsgw_vxc(
                driver_params.input_dir, driver_params.prefix_vxc_band,
                *band_reference, aims_input, true, basis);
    });

    HamiltonianCutOptions cut_options;
    cut_options.unoccupied_keep = driver_params.qsgw_unoccupied_keep;
    cut_options.mode = hamiltonian_cut_mode_from_int(driver_params.qsgw_cut_mode);
    cut_options.shift_ha = driver_params.qsgw_cut_shift_ha;
    const bool apply_cut = cut_options.mode != HamiltonianCutMode::Uncut;
    SpinKMatrixMap current_hamiltonian =
        apply_cut ? apply_hamiltonian_cut(reference_hamiltonian, reference_hamiltonian, dataset->mf,
                                          cut_options)
                  : reference_hamiltonian;
    SpinKMatrixMap current_band_hamiltonian = compute_band
        ? apply_hamiltonian_cut(
              band_reference_hamiltonian, band_reference_hamiltonian,
              dataset->mf_band, cut_options)
        : band_reference_hamiltonian;

    const CheckpointContext checkpoint_context{reference,
                                               band_reference ? &*band_reference : nullptr,
                                               dataset->pbc.kfrac_list,
                                               dataset->kfrac_band_list,
                                               dataset->pbc.weight_k,
                                               dft_vxc,
                                               dft_vxc_band,
                                               driver_params.qsgw_cut_mode,
                                               driver_params.qsgw_unoccupied_keep,
                                               driver_params.qsgw_cut_shift_ha,
                                               mixing_beta,
                                               opts.nfreq,
                                               update_head ? opts.option_dielect_func : 0,
                                               driver_params.qsgw_headwing_velocity,
                                               ao_operators ? &*ao_operators : nullptr,
                                               &reference_velocity,
                                               &dataset->pbc.period,
                                               &dataset->pbc.latvec};
    int starting_iteration = 0;
    IterationSummary starting_summary;
    starting_summary.fermi_energy_ev = reference.get_efermi() * HA2EV;
    starting_summary.gap_ev = initial_occupations.gap * HA2EV;
    starting_summary.electron_count = initial_occupations.electron_count;
    starting_summary.beta = mixing_beta;
    if (driver_params.qsgw_restart_file != "none")
    {
        collective_root_stage(
            dataset->comm_h, "QSGW checkpoint restore",
            [&]
            {
                starting_summary =
                    read_qsgw_checkpoint(driver_params.qsgw_restart_file, checkpoint_context,
                                         dataset->mf, compute_band ? &dataset->mf_band : nullptr,
                                         current_hamiltonian, current_band_hamiltonian);
                starting_iteration = starting_summary.iteration;
                if (starting_iteration >= driver_params.qsgw_max_iter)
                    throw std::runtime_error(
                        "qsgw_max_iter must exceed the checkpoint's completed iteration");
                lib_printf("QSGW restored completed iteration %d from %s\n", starting_iteration,
                           driver_params.qsgw_restart_file.c_str());
            });
        dataset->comm_h.bcast(&starting_iteration, 1, 0);
        broadcast_checkpoint_meanfield(dataset->mf, dataset->comm_h);
        broadcast_spin_k_matrix_map(current_hamiltonian, 0, dataset->comm_h);
        if (compute_band)
        {
            broadcast_checkpoint_meanfield(dataset->mf_band, dataset->comm_h);
            broadcast_spin_k_matrix_map(current_band_hamiltonian, 0, dataset->comm_h);
            dataset->is_band_calc_done = true;
        }
    }

    std::ofstream qp_output;
    SpinKMatrixMap qp_rotation, qp_band_rotation;
    if (driver_params.qsgw_output_qp)
        collective_root_stage(dataset->comm_h, "QSGW QP frame recovery",
                              [&]
                              {
                                  qp_rotation = reference_basis_rotation(dataset->mf, reference);
                                  if (compute_band)
                                      qp_band_rotation = reference_basis_rotation(dataset->mf_band,
                                                                                  *band_reference);
                              });

    std::ofstream trace;
    std::ofstream eigenvalue_trace;
    std::ofstream matrix_trace;
    collective_root_stage(
        dataset->comm_h, "QSGW trace setup",
        [&]
        {
            if (driver_params.qsgw_output_qp)
            {
                qp_output.open(join_path(opts.output_dir, "qsgw_qp.csv"));
                if (!qp_output) throw LIBRPA_RUNTIME_ERROR("Cannot open QSGW quasiparticle output");
                qp_output << "state_index,gw_update,channel,spin,kpoint,kx,ky,kz,band,occupation,"
                             "e0_eV,vstatic_eV,sigma_x_eV,re_sigma_c_effective_eV,im_sigma_c_eV,"
                             "eqp_eV,solver_info\n"
                          << std::scientific << std::setprecision(17);
            }
            trace.open(join_path(opts.output_dir, "qsgw_iterations.dat"));
            eigenvalue_trace.open(join_path(opts.output_dir, "qsgw_eigenvalues.dat"));
            if (driver_params.qsgw_write_iteration_matrices)
                matrix_trace.open(join_path(opts.output_dir, "qsgw_matrices.dat"));
            if (!trace || !eigenvalue_trace ||
                (driver_params.qsgw_write_iteration_matrices && !matrix_trace))
                throw std::runtime_error("Cannot open QSGW trace output");
            write_iteration_summary_header(trace);
            write_eigenvalue_trace_header(eigenvalue_trace);
            if (driver_params.qsgw_write_iteration_matrices)
            {
                write_matrix_trace_header(matrix_trace);
                write_matrix_component_trace(matrix_trace, 0, IterationChannel::Grid, "h0",
                                             reference_hamiltonian);
                write_matrix_component_trace(matrix_trace, 0, IterationChannel::Grid, "vxc_dft",
                                             dft_vxc);
                write_wavefunction_trace(matrix_trace, 0, IterationChannel::Grid, reference);
                write_occupation_trace(matrix_trace, 0, IterationChannel::Grid, reference);
                if (compute_band)
                {
                    write_matrix_component_trace(matrix_trace, 0, IterationChannel::Band, "h0",
                                                 band_reference_hamiltonian);
                    write_matrix_component_trace(matrix_trace, 0, IterationChannel::Band, "vxc_dft",
                                                 dft_vxc_band);
                    write_wavefunction_trace(matrix_trace, 0, IterationChannel::Band,
                                             *band_reference);
                    write_occupation_trace(matrix_trace, 0, IterationChannel::Band,
                                           *band_reference);
                }
                if (starting_iteration > 0)
                {
                    write_matrix_component_trace(matrix_trace, starting_iteration,
                                                 IterationChannel::Grid, "mixed_h",
                                                 current_hamiltonian);
                    write_wavefunction_trace(matrix_trace, starting_iteration,
                                             IterationChannel::Grid, dataset->mf);
                    write_occupation_trace(matrix_trace, starting_iteration, IterationChannel::Grid,
                                           dataset->mf);
                    if (compute_band)
                    {
                        write_matrix_component_trace(matrix_trace, starting_iteration,
                                                     IterationChannel::Band, "mixed_h",
                                                     current_band_hamiltonian);
                        write_wavefunction_trace(matrix_trace, starting_iteration,
                                                 IterationChannel::Band, dataset->mf_band);
                        write_occupation_trace(matrix_trace, starting_iteration,
                                               IterationChannel::Band, dataset->mf_band);
                    }
                }
            }
            write_iteration_summary(trace, starting_summary);
            write_eigenvalue_trace(eigenvalue_trace, starting_iteration, IterationChannel::Grid,
                                   dataset->mf, dataset->pbc.kfrac_list);
            if (compute_band)
            {
                write_eigenvalue_trace(eigenvalue_trace, starting_iteration, IterationChannel::Band,
                                       dataset->mf_band, dataset->kfrac_band_list);
            }
        });

    bool converged = false;
    int completed_iterations = starting_iteration;
    const auto export_hamiltonian = [&](const int iteration, const SpinKMatrixMap& hamiltonian)
    {
        if (!driver_params.qsgw_export_hamiltonian_for_pyatb) return;
        collective_root_stage(
            dataset->comm_h, "QSGW complete Hamiltonian HR export",
            [&]
            {
                const auto exported = reconstruct_qsgw_hamiltonian(
                    reference, hamiltonian, *ao_operators, dataset->pbc.kfrac_list,
                    dataset->pbc.period, dataset->pbc.latvec);
                for (int spin = 0; spin < reference.get_n_spins(); ++spin)
                {
                    std::ostringstream name;
                    name << "hrs" << spin + 1 << "_nao_qsgw_iter_" << std::setfill('0')
                         << std::setw(4) << iteration << ".csr";
                    driver::qsgw::write_qsgw_hamiltonian_csr(
                        join_path(opts.output_dir, name.str()), exported.hamiltonian.at(spin),
                        *ao_operators, iteration, reference.get_n_spins(), spin + 1);
                }
                lib_printf("QSGW HR iteration %d: full complex Ry, metric %.8e FT_Ha %.8e\n",
                           iteration, exported.reference_metric_residual,
                           exported.fourier_roundtrip_residual_ha);
            });
    };
    // Validate export inputs before the first GW build and retain the actual
    // entering state. Fresh iteration 0 is PBE; a restart keeps its saved index.
    export_hamiltonian(starting_iteration,
                       starting_iteration == 0 ? reference_hamiltonian : current_hamiltonian);
    for (int iteration = starting_iteration + 1; iteration <= driver_params.qsgw_max_iter;
         ++iteration)
    {
        const int entering_state = iteration - 1;
        completed_iterations = iteration;
        const EigenvalueSnapshot previous = eigenvalue_snapshot(dataset->mf);
        dataset->invalidate_compute_objects();
        if (update_head)
        {
            collective_all_stage(
                dataset->comm_h, "QSGW live head/wing velocity",
                [&]
                {
                    if (driver_params.qsgw_headwing_velocity == "current_hr")
                    {
                        auto velocity = compute_live_headwing_velocity(
                            reference, dataset->mf,
                            // The initial cut H is a residual reference; update 1
                            // still screens with the actual entering PBE eigenpairs.
                            starting_iteration == 0 && iteration == 1 ? reference_hamiltonian
                                                                      : current_hamiltonian,
                            *ao_operators, dataset->pbc.kfrac_list, dataset->pbc.period,
                            dataset->pbc.latvec);
                        dataset->velocity_matrix = std::move(velocity.velocity);
                        if (dataset->comm_h.is_root())
                            lib_printf(
                                "QSGW update %d current-HR velocity: metric %.8e live_metric %.8e "
                                "frame_Ha %.8e FT_Ha %.8e\n",
                                iteration, velocity.reference_metric_residual,
                                velocity.live_metric_residual,
                                velocity.hamiltonian_frame_residual_ha,
                                velocity.fourier_roundtrip_residual_ha);
                    }
                    else
                    {
                        dataset->velocity_matrix = reference_velocity;
                        if (driver_params.qsgw_headwing_velocity == "rotate_reference")
                            align_velocity_to_reference_wfc(reference, dataset->mf,
                                                            dataset->velocity_matrix);
                    }
                });
            dataset->p_headwing.reset();
            initialize_ds_tfgrids(*dataset, opts);
            // Legacy option-4 head reused the first k-point occupation row.
            // Preserve it only for the unchanged reference-entries map.
            if (driver_params.qsgw_headwing_velocity == "reference")
            {
                ScopedLegacyHeadOccupations legacy_head_occupations(dataset->mf);
                initialize_ds_headwing(*dataset, opts, false);
            }
            else
            {
                initialize_ds_headwing(*dataset, opts, update_wing);
            }
            // The analytic object is replaced every update. Reapply the fixed producer
            // normalization while retaining its live eigenpairs, velocity and frequencies.
            collective_all_stage(dataset->comm_h, "QSGW head/wing Coulomb normalization",
                                 [&]
                                 {
                                     configure_headwing_coulomb_head(driver_params.input_dir,
                                                                     *dataset->p_headwing);
                                 });
            if (update_wing && (!dataset->p_headwing || !dataset->p_headwing->has_wing()))
                throw LIBRPA_RUNTIME_ERROR("QSGW analytic wing was not rebuilt for the live state");
            collective_root_stage(
                dataset->comm_h, "QSGW head/wing numerical trace",
                [&]
                {
                    if (!driver_params.qsgw_write_iteration_matrices ||
                        driver_params.qsgw_headwing_velocity == "reference")
                        return;
                    for (int a = 0; a < 3; ++a)
                    {
                        SpinKMatrixMap velocities;
                        for (int spin = 0; spin < dataset->mf.get_n_spins(); ++spin)
                            for (int ik = 0; ik < dataset->mf.get_n_kpoints(); ++ik)
                            {
                                const auto& v = dataset->velocity_matrix[spin][ik][a];
                                Matz matrix(v.nr, v.nc, MAJOR::ROW);
                                for (int i = 0; i < v.nr; ++i)
                                    for (int j = 0; j < v.nc; ++j) matrix(i, j) = v(i, j);
                                velocities[spin][ik] = std::move(matrix);
                            }
                        write_matrix_component_trace(
                            matrix_trace, iteration, IterationChannel::Grid,
                            "head_velocity_" + std::to_string(a), velocities);
                    }
                    const auto& frequencies = dataset->p_headwing->get_frequency_nodes();
                    for (const auto& [name, tensors] :
                         std::vector<std::pair<std::string, const std::vector<Matz>*>>{
                             {"head_tensor", &dataset->p_headwing->get_head_tensors()},
                             {"wing_auxiliary", &dataset->p_headwing->get_wing_auxiliary()}})
                    {
                        SpinKFrequencyMatrixMap matrices;
                        for (std::size_t i = 0; i < tensors->size(); ++i)
                            matrices[0][0][frequencies.at(i)] = tensors->at(i);
                        write_frequency_matrix_component_trace(
                            matrix_trace, iteration, IterationChannel::Grid, name, matrices);
                    }
                });
        }
        h.build_g0w0_sigma(opts);
        if (!dataset->p_exx || !dataset->p_g0w0)
            throw LIBRPA_RUNTIME_ERROR(
                "QSGW failed to build live EXX/Sigma real-space objects");
        if (update_head && !dataset->p_headwing)
            throw LIBRPA_RUNTIME_ERROR(
                "QSGW head-only object was not rebuilt on the live frequency grid");

        {
            // The upstream projection path uses the output flag as its
            // full-matrix retention gate. Restore it before leaving this
            // scope so QSGW does not change the user-facing output setting.
            ScopedSigmaMatrixRetention retain_sigma_matrices(
                dataset->p_g0w0->output_sigc_ks_mat_kf);
            ScopedReferenceEigenvectors fixed_basis_projection(
                dataset->mf, reference);
            dataset->p_exx->build_KS_kgrid_blacs(
                dataset->blacs_h,
                opts.use_gpu_replace_scalapack == LIBRPA_SWITCH_ON);
            dataset->p_g0w0->build_sigc_matrix_KS_kgrid_blacs(
                dataset->blacs_h,
                opts.use_gpu_replace_scalapack == LIBRPA_SWITCH_ON);
        }
        const std::vector<double> frequencies =
            dataset->tfg.get_freq_nodes();
        const SigmaMatrixMap sigma = collect_sigma_root(
            *dataset->p_g0w0, reference, frequencies, dataset->comm_h);

        SpinKMatrixMap exchange;
        collective_root_stage(dataset->comm_h, "QSGW EXX collection", [&] {
            exchange = copy_exx_root(*dataset->p_exx, reference);
        });

        SigmaMatrixMap sigma_band;
        SpinKMatrixMap exchange_band;
        if (compute_band)
        {
            dataset->p_exx->reset_kspace();
            dataset->p_g0w0->reset_kspace();
            const auto bvk_remap =
                librpa_int::api::build_band_bvk_remap(
                    dataset->atoms, dataset->pbc, opts.option_bvk_remap);
            {
                ScopedSigmaMatrixRetention retain_sigma_matrices(
                    dataset->p_g0w0->output_sigc_ks_mat_kf);
                dataset->p_exx->build_KS_band_blacs(
                    band_reference->get_eigenvectors(),
                    dataset->kfrac_band_list, bvk_remap,
                    dataset->blacs_h,
                    opts.use_gpu_replace_scalapack ==
                        LIBRPA_SWITCH_ON);
                dataset->p_g0w0->build_sigc_matrix_KS_band_blacs(
                    band_reference->get_eigenvectors(),
                    dataset->kfrac_band_list, bvk_remap,
                    dataset->blacs_h,
                    opts.use_gpu_replace_scalapack ==
                        LIBRPA_SWITCH_ON);
            }
            sigma_band = collect_sigma_root(
                *dataset->p_g0w0, *band_reference,
                frequencies, dataset->comm_h);
            collective_root_stage(
                dataset->comm_h, "QSGW band EXX collection", [&] {
                    exchange_band = copy_exx_root(
                        *dataset->p_exx, *band_reference);
                });
        }

        if (driver_params.qsgw_output_qp)
            collective_root_stage(
                dataset->comm_h, "QSGW entering-state QP equation",
                [&]
                {
                    write_qsgw_qp_state(
                        qp_output, entering_state, iteration, "grid", dataset->mf,
                        dataset->pbc.kfrac_list, qp_rotation,
                        entering_state == 0 ? reference_hamiltonian : current_hamiltonian,
                        reference_hamiltonian, dft_vxc, exchange, sigma, dataset->tfg,
                        dataset->pbc.weight_k);
                    if (compute_band)
                        write_qsgw_qp_state(qp_output, entering_state, iteration, "band",
                                            dataset->mf_band, dataset->kfrac_band_list,
                                            qp_band_rotation,
                                            entering_state == 0 ? band_reference_hamiltonian
                                                                : current_band_hamiltonian,
                                            band_reference_hamiltonian, dft_vxc_band, exchange_band,
                                            sigma_band, dataset->tfg, {});
                });

        SpinKMatrixMap mixed_hamiltonian;
        SpinKMatrixMap mixed_band_hamiltonian;
        SpinKMatrixMap band_exchange_output;
        double residual_l2 = 0.0;
        double residual_max = 0.0;
        std::string matrix_rows;
        collective_root_stage(dataset->comm_h, "QSGW Hamiltonian update", [&] {
            const SpinKMatrixMap correlation = build_correlation_map(
                dataset->mf, sigma, frequencies, opts);
            SpinKMatrixMap correlation_band;
            if (compute_band)
            {
                correlation_band = build_correlation_map(
                    dataset->mf_band, sigma_band, frequencies, opts);
            }
            const SpinKMatrixMap raw_uncut = assemble_effective_hamiltonian(
                reference_hamiltonian, dft_vxc, exchange, correlation);
            const SpinKMatrixMap raw = apply_cut
                                           ? apply_hamiltonian_cut(raw_uncut, reference_hamiltonian,
                                                                   dataset->mf, cut_options)
                                           : raw_uncut;
            SpinKMatrixMap raw_band;
            if (compute_band)
            {
                band_exchange_output = exchange_band;
                const SpinKMatrixMap raw_band_uncut =
                    assemble_effective_hamiltonian(
                        band_reference_hamiltonian, dft_vxc_band,
                        exchange_band, correlation_band);
                raw_band = apply_hamiltonian_cut(
                    raw_band_uncut, band_reference_hamiltonian,
                    dataset->mf_band, cut_options);
            }
            const auto residual = measure_spin_k_hamiltonian_residual(
                raw, current_hamiltonian);
            residual_l2 = residual.l2;
            residual_max = residual.maximum;
            mixed_hamiltonian = mix_spin_k_hamiltonian(
                current_hamiltonian, raw, mixing_beta);
            if (compute_band)
                mixed_band_hamiltonian = mix_spin_k_hamiltonian(
                    current_band_hamiltonian, raw_band, mixing_beta);
            if (apply_cut)
            {
                mixed_hamiltonian = apply_hamiltonian_cut(
                    mixed_hamiltonian, reference_hamiltonian, dataset->mf,
                    cut_options);
                if (compute_band)
                    mixed_band_hamiltonian =
                        apply_hamiltonian_cut(mixed_band_hamiltonian, band_reference_hamiltonian,
                                              dataset->mf_band, cut_options);
            }
            if (driver_params.qsgw_write_iteration_matrices)
            {
                std::ostringstream rows;
                write_frequency_matrix_component_trace(
                    rows, iteration, IterationChannel::Grid,
                    "sigma_c_iw", sigma);
                write_matrix_component_trace(
                    rows, iteration, IterationChannel::Grid, "exx",
                    exchange);
                write_matrix_component_trace(
                    rows, iteration, IterationChannel::Grid, "vc",
                    correlation);
                write_matrix_component_trace(
                    rows, iteration, IterationChannel::Grid, "raw_h", raw);
                write_matrix_component_trace(
                    rows, iteration, IterationChannel::Grid, "mixed_h",
                    mixed_hamiltonian);
                if (compute_band)
                {
                    write_matrix_component_trace(
                        rows, iteration, IterationChannel::Band, "exx",
                        exchange_band);
                    write_matrix_component_trace(
                        rows, iteration, IterationChannel::Band, "vc",
                        correlation_band);
                    write_matrix_component_trace(
                        rows, iteration, IterationChannel::Band, "raw_h",
                        raw_band);
                    write_matrix_component_trace(
                        rows, iteration, IterationChannel::Band, "mixed_h",
                        mixed_band_hamiltonian);
                }
                matrix_rows = rows.str();
            }
        });

        broadcast_spin_k_matrix_map(
            mixed_hamiltonian, 0, dataset->comm_h);
        if (compute_band)
        {
            broadcast_spin_k_matrix_map(
                mixed_band_hamiltonian, 0, dataset->comm_h);
        }
        // The next velocity is evaluated on every rank. Its Hamiltonian must
        // advance with the live eigenpairs, not remain a root-only residual map.
        current_hamiltonian = mixed_hamiltonian;
        if (compute_band) current_band_hamiltonian = mixed_band_hamiltonian;
        // Export the exact mixed/post-cut operator about to be diagonalized.
        // AO reconstruction uses the immutable initial coefficients, so the
        // eigenpairs produced below are not needed for this output.
        export_hamiltonian(iteration, current_hamiltonian);
        const FixedBasisDiagonalizationResult diagonalization =
            diagonalize_in_reference_basis(
                dataset->mf, reference, mixed_hamiltonian);
        std::optional<FixedBasisDiagonalizationResult> band_diagonalization;
        if (compute_band)
        {
            band_diagonalization = diagonalize_in_reference_basis(
                dataset->mf_band, *band_reference,
                mixed_band_hamiltonian);
            dataset->is_band_calc_done = true;
        }
        const OccupationResult occupations = update_qsgw_occupations(
            dataset->mf, reference, dataset->pbc.weight_k, electron_count);
        if (driver_params.qsgw_output_qp && dataset->comm_h.is_root())
        {
            qp_rotation = diagonalization.unitary;
            if (compute_band) qp_band_rotation = band_diagonalization->unitary;
        }
        if (compute_band)
            dataset->mf_band.get_efermi() = occupations.chemical_potential;
        const double maximum_change_ev =
            max_eigenvalue_change(dataset->mf, previous) * HA2EV;

        int converged_flag = 0;
        collective_root_stage(dataset->comm_h, "QSGW reporting", [&] {
            converged_flag = qsgw_iteration_converged(
                                 iteration, driver_params.qsgw_min_iter,
                                 maximum_change_ev,
                                 driver_params.qsgw_convergence_tolerance_ev)
                                 ? 1
                                 : 0;
            IterationSummary summary;
            summary.iteration = iteration;
            summary.maximum_eigenvalue_change_ev = maximum_change_ev;
            summary.residual_l2_ha = residual_l2;
            summary.residual_max_ha = residual_max;
            summary.fermi_energy_ev =
                occupations.chemical_potential * HA2EV;
            summary.gap_ev = occupations.gap * HA2EV;
            summary.electron_count = occupations.electron_count;
            summary.converged = converged_flag != 0;
            summary.beta = mixing_beta;
            if (driver_params.qsgw_checkpoint_file != "none")
                write_qsgw_checkpoint(
                    resolve_input_path(opts.output_dir, driver_params.qsgw_checkpoint_file),
                    checkpoint_context, summary, dataset->mf,
                    compute_band ? &dataset->mf_band : nullptr, current_hamiltonian,
                    current_band_hamiltonian);
            write_iteration_summary(trace, summary);
            write_eigenvalue_trace(
                eigenvalue_trace, iteration, IterationChannel::Grid,
                dataset->mf, dataset->pbc.kfrac_list);
            if (compute_band)
            {
                write_eigenvalue_trace(
                    eigenvalue_trace, iteration, IterationChannel::Band,
                    dataset->mf_band, dataset->kfrac_band_list);
                for (int spin = 0;
                     spin < dataset->mf_band.get_n_spins(); ++spin)
                {
                    const auto make_filename = [&](const std::string& prefix) {
                        std::ostringstream name;
                        name << prefix << spin + 1 << '_' << iteration
                             << ".dat";
                        return join_path(opts.output_dir, name.str());
                    };
                    std::ofstream ks_output(
                        make_filename("KS_band_spin_"));
                    std::ofstream exx_output(
                        make_filename("EXX_band_spin_"));
                    std::ofstream qsgw_output(
                        make_filename("QSGW_band_spin_"));
                    write_qsgw_band_spin_tables(
                        ks_output, exx_output, qsgw_output,
                        dataset->mf_band, *band_reference,
                        dataset->kfrac_band_list,
                        band_reference_hamiltonian, dft_vxc_band,
                        band_exchange_output, spin,
                        occupations.chemical_potential);
                }
            }
            if (driver_params.qsgw_write_iteration_matrices)
            {
                matrix_trace << matrix_rows;
                write_matrix_component_trace(
                    matrix_trace, iteration, IterationChannel::Grid,
                    "rotation_u", diagonalization.unitary);
                write_wavefunction_trace(
                    matrix_trace, iteration, IterationChannel::Grid,
                    dataset->mf);
                write_occupation_trace(
                    matrix_trace, iteration, IterationChannel::Grid,
                    dataset->mf);
                if (compute_band)
                {
                    write_matrix_component_trace(
                        matrix_trace, iteration, IterationChannel::Band,
                        "rotation_u", band_diagonalization->unitary);
                    write_wavefunction_trace(
                        matrix_trace, iteration, IterationChannel::Band,
                        dataset->mf_band);
                    write_occupation_trace(
                        matrix_trace, iteration, IterationChannel::Band,
                        dataset->mf_band);
                }
            }
            trace.flush();
            eigenvalue_trace.flush();
            if (driver_params.qsgw_write_iteration_matrices)
                matrix_trace.flush();
            lib_printf("QSGW iteration %d: max_delta=% .8e eV gap=% .8f eV mixer=%s%s\n", iteration,
                       maximum_change_ev, occupations.gap * HA2EV, effective_mixer.c_str(),
                       converged_flag ? ", converged" : "");
        });
        dataset->comm_h.bcast(&converged_flag, 1, 0);
        converged = converged_flag != 0;
        if (converged) break;
    }

    if (dataset->comm_h.is_root())
    {
        if (!converged)
            lib_printf(LIBRPA_VERBOSE_WARN,
                       "QSGW reached qsgw_max_iter=%d without convergence\n",
                       driver_params.qsgw_max_iter);
        lib_printf("QSGW completed iterations: %d\n", completed_iterations);
    }
    profiler.stop("qsgw");
}

} // namespace

void driver::task_qsgw() { run_qsgw(); }
