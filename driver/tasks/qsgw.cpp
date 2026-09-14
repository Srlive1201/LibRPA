#include <algorithm>
#include <cmath>
#include <exception>
#include <fstream>
#include <iomanip>
#include <map>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <librpa_enums.h>

#include "../../src/api/compute_helper.h"
#include "../../src/api/dataset_helper.h"
#include "../../src/api/instance_manager.h"
#include "../../src/io/fs.h"
#include "../../src/io/global_io.h"
#include "../qsgw/band_output.h"
#include "../../src/qsgw/convergence.h"
#include "../../src/qsgw/correlation_potential.h"
#include "../../src/qsgw/distributed_matrix.h"
#include "../../src/qsgw/effective_hamiltonian.h"
#include "../../src/qsgw/fixed_basis.h"
#include "../../src/qsgw/hamiltonian_cut.h"
#include "../../src/qsgw/hamiltonian_mixing.h"
#include "../qsgw/iteration_trace.h"
#include "../../src/qsgw/occupation.h"
#include "../../src/qsgw/projection_target.h"
#include "../qsgw/vxc_io.h"
#include "../../src/utils/constants.h"
#include "../../src/utils/profiler.h"
#include "../driver.h"
#include "../read_data.h"
#include "../reader_coulomb.h"
#include "../task.h"

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

void prepare_stage_one_symmetry_context(librpa_int::Dataset& dataset)
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
    using librpa_int::qsgw::collect_blacs_matrix_root;
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

void run_qsgw_stage_one(const bool compute_band)
{
    using namespace driver;
    using namespace librpa_int;
    using namespace librpa_int::global;
    using namespace librpa_int::qsgw;

    const bool update_head =
        driver::get_bool(opts.replace_w_head) &&
        opts.option_dielect_func == 4;
    if (driver_params.use_pyatb)
    {
        throw LIBRPA_RUNTIME_ERROR(
            "QSGW independent PyATB head updates are unsupported; use the same-grid velocity input");
    }
    if (driver::get_bool(opts.replace_w_head) && !update_head)
    {
        throw LIBRPA_RUNTIME_ERROR(
            "QSGW supports only analytic head-only mode; set option_dielect_func = 4");
    }
    if (driver::get_bool(opts.use_kpara_scf_eigvec))
    {
        throw LIBRPA_RUNTIME_ERROR(
            "QSGW fixed-basis iteration requires replicated SCF wavefunctions; set use_kpara_scf_eigvec = false");
    }

    profiler.start("qsgw", "QSGW fixed-basis self-consistent calculation");
    const auto dataset = api::get_dataset_instance(h);
    if (opts.parallel_routing != LIBRPA_ROUTING_LIBRI &&
        opts.parallel_routing != LIBRPA_ROUTING_AUTO)
    {
        throw LIBRPA_RUNTIME_ERROR(
            "QSGW stage one requires parallel_routing=libri or auto");
    }
    if (compute_band)
    {
        const std::string band_kpath_path = resolve_input_path(
            driver_params.input_dir, driver_params.fn_band_kpath_info);
        read_band_kpath_info(band_kpath_path);
        dataset->comm_h.barrier();
        read_band_meanfield_data(driver_params.input_dir);
        dataset->comm_h.barrier();
    }
    read_Vq_row(driver_params.input_dir, driver_params.prefix_coul_cut,
                opts.vq_threshold, local_atpair, true,
                driver_params.version_coul_reader,
                driver::get_bool(opts.use_shrink_abfs));
    collective_root_stage(dataset->comm_h, "QSGW reference validation", [&] {
        prepare_stage_one_symmetry_context(*dataset);
        validate_qsgw_reference(*dataset, compute_band);
    });

    const MeanField reference = dataset->mf;
    const double electron_count = physical_electron_count(
        reference, dataset->pbc.weight_k);
    const OccupationResult initial_occupations = analyze_qsgw_occupations(
        reference, dataset->pbc.weight_k, electron_count);
    const SpinKMatrixMap reference_hamiltonian =
        build_reference_hamiltonian(reference);

    VelocityMatrix reference_velocity;
    if (update_head)
    {
        read_headwing_input(driver_params.input_dir, false);
        reference_velocity = dataset->velocity_matrix;
        if (driver_params.constants_choice == "aims")
        {
            prepare_fhi_aims_interband_velocity(
                reference_velocity, reference);
            align_distributed_velocity_to_reference_wfc(
                dataset->p_headwing->get_meanfield_df(), reference,
                reference_velocity, dataset->comm_h);
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
    cut_options.unoccupied_keep =
        driver_params.qsgw_band0_unoccupied_keep;
    cut_options.mode = hamiltonian_cut_mode_from_int(
        driver_params.qsgw_band0_cut_mode);
    cut_options.shift_ha = driver_params.qsgw_band0_cut_shift_ha;
    SpinKMatrixMap current_hamiltonian = compute_band
        ? apply_hamiltonian_cut(
              reference_hamiltonian, reference_hamiltonian, dataset->mf,
              cut_options)
        : reference_hamiltonian;
    SpinKMatrixMap current_band_hamiltonian = compute_band
        ? apply_hamiltonian_cut(
              band_reference_hamiltonian, band_reference_hamiltonian,
              dataset->mf_band, cut_options)
        : band_reference_hamiltonian;
    const double mixing_beta = driver_params.qsgw_mixer == "linear"
        ? driver_params.qsgw_mixing_beta : 1.0;

    std::ofstream trace;
    std::ofstream eigenvalue_trace;
    std::ofstream matrix_trace;
    collective_root_stage(dataset->comm_h, "QSGW trace setup", [&] {
        trace.open(join_path(opts.output_dir, "qsgw_iterations.dat"));
        eigenvalue_trace.open(
            join_path(opts.output_dir, "qsgw_eigenvalues.dat"));
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
            write_matrix_component_trace(
                matrix_trace, 0, IterationChannel::Grid, "h0",
                reference_hamiltonian);
            write_matrix_component_trace(
                matrix_trace, 0, IterationChannel::Grid, "vxc_dft",
                dft_vxc);
            write_wavefunction_trace(
                matrix_trace, 0, IterationChannel::Grid, reference);
            write_occupation_trace(
                matrix_trace, 0, IterationChannel::Grid, dataset->mf);
            if (compute_band)
            {
                write_matrix_component_trace(
                    matrix_trace, 0, IterationChannel::Band, "h0",
                    band_reference_hamiltonian);
                write_matrix_component_trace(
                    matrix_trace, 0, IterationChannel::Band, "vxc_dft",
                    dft_vxc_band);
                write_wavefunction_trace(
                    matrix_trace, 0, IterationChannel::Band,
                    *band_reference);
                write_occupation_trace(
                    matrix_trace, 0, IterationChannel::Band,
                    dataset->mf_band);
            }
        }
        IterationSummary summary;
        summary.iteration = 0;
        summary.fermi_energy_ev = reference.get_efermi() * HA2EV;
        summary.gap_ev = initial_occupations.gap * HA2EV;
        summary.electron_count = initial_occupations.electron_count;
        summary.beta = mixing_beta;
        write_iteration_summary(trace, summary);
        write_eigenvalue_trace(eigenvalue_trace, 0,
                               IterationChannel::Grid, reference,
                               dataset->pbc.kfrac_list);
        if (compute_band)
        {
            write_eigenvalue_trace(
                eigenvalue_trace, 0, IterationChannel::Band,
                *band_reference, dataset->kfrac_band_list);
        }
    });

    bool converged = false;
    int completed_iterations = 0;
    for (int iteration = 1;
         iteration <= driver_params.qsgw_max_iter; ++iteration)
    {
        completed_iterations = iteration;
        const EigenvalueSnapshot previous = eigenvalue_snapshot(dataset->mf);
        dataset->invalidate_compute_objects();
        if (update_head)
        {
            dataset->velocity_matrix = reference_velocity;
            dataset->p_headwing.reset();
            initialize_ds_tfgrids(*dataset, opts);
            // Legacy option-4 head reused the first k-point occupation row.
            // Limit that compatibility behavior to head construction.
            ScopedLegacyHeadOccupations legacy_head_occupations(dataset->mf);
            initialize_ds_headwing(*dataset, opts, false);
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
            const SpinKMatrixMap raw = compute_band
                ? apply_hamiltonian_cut(
                      raw_uncut, reference_hamiltonian, dataset->mf,
                      cut_options)
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
            if (compute_band)
            {
                mixed_hamiltonian = apply_hamiltonian_cut(
                    mixed_hamiltonian, reference_hamiltonian, dataset->mf,
                    cut_options);
                mixed_band_hamiltonian = apply_hamiltonian_cut(
                    mixed_band_hamiltonian, band_reference_hamiltonian,
                    dataset->mf_band, cut_options);
            }
            current_hamiltonian = mixed_hamiltonian;
            if (compute_band)
                current_band_hamiltonian = mixed_band_hamiltonian;

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
            lib_printf(
                "QSGW iteration %d: max_delta=% .8e eV gap=% .8f eV mixer=%s%s\n",
                iteration, maximum_change_ev, occupations.gap * HA2EV,
                driver_params.qsgw_mixer.c_str(),
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

void driver::task_qsgw()
{
    run_qsgw_stage_one(false);
}

void driver::task_qsgw_band()
{
    run_qsgw_stage_one(true);
}
