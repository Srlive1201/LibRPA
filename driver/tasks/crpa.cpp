#include "crpa.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>

#include "../../src/io/fs.h"
#include "../../src/io/global_io.h"
#include "../../src/mpi/global_mpi.h"
#include "../../src/utils/constants.h"
#include "../driver.h"
#include "../read_data.h"
#include "../reader/reader_nao_overlap.h"
#include "../task.h"
#include "crpa_model_output.h"
#include "librpa_crpa.h"
#include "librpa_input.h"
#include "librpa_orbital.h"

namespace
{

std::ofstream output_file(const std::string& directory, const std::string& name)
{
    std::ofstream output(librpa_int::join_path(directory, name));
    if (!output) throw std::runtime_error("cannot create cRPA output " + name);
    output.exceptions(std::ios::badbit | std::ios::failbit);
    output << std::setprecision(16);
    return output;
}

void write_orbital_analysis(const LibrpaOrbitalResult* result,
                            const std::vector<std::string>& species, const std::string& directory)
{
    using librpa_int::HA2EV;
    const auto& params = driver::driver_params;
    LibrpaOrbitalInfo info;
    librpa_orbital_info(result, &info);
    const double reference =
        params.crpa_pdos_reference == "fermi" ? info.fermi_ha : info.highest_occupied_ha;
    if (!std::isfinite(reference))
        throw std::runtime_error(
            "no finite PDOS energy reference; use fermi if no occupied state exists");
    auto metadata =
        output_file(directory, params.crpa_output_format == "dat" ? "crpa_orbital_analysis.dat"
                                                                  : "crpa_orbital_analysis.txt");
    metadata
        << "Existing scalar/spinor KS states; no Hamiltonian read or physical KS solve.\n"
        << "Mulliken=sum_s Re[c_mu,s* (S c_s)_mu]; "
           "Lowdin=sum_s |(S^(1/2)c_s)_mu|^2/(sum_s c_s-dagger S c_s).\n"
        << (info.n_spinor == 2 ? "Spinor components are summed into one stored spin channel; "
                                 "its index 0 is not a conserved spin projection.\n"
                               : "Scalar inputs retain their separate stored spin channels.\n")
        << "Groups include every radial shell; eg/t2g use global Cartesian axes.\n"
        << "PDOS uses loaded k weights and separate stored spin channels; no occupation or extra "
           "degeneracy factor.\n"
        << "PDOS character is diagnostic; orbital_weights uses the actual local output projector.\n"
        << "energy_reference = " << params.crpa_pdos_reference << '\n'
        << "energy_reference_absolute_ev = " << reference * HA2EV << '\n'
        << "sigma_ev = " << params.crpa_pdos_sigma_ev << '\n'
        << "n_spins = " << info.n_spins << '\n'
        << "n_spinor = " << info.n_spinor << '\n'
        << "n_kpoints = " << info.n_kpoints << '\n'
        << "n_bands = " << info.n_bands << '\n'
        << "n_aos = " << info.n_aos << '\n'
        << "min_overlap_eigenvalue = " << info.min_overlap_eigenvalue << '\n'
        << "max_KS_norm_error = " << info.max_norm_error << '\n'
        << "max_Lowdin_sum_error = " << info.max_sum_error << '\n'
        << "max_overlap_sqrt_residual = " << info.max_sqrt_residual << '\n';

    auto groups = output_file(directory, "crpa_orbital_groups.dat");
    groups << "# group atom species l channel; indices zero based, all radial shells summed\n";
    for (int g = 0; g < info.n_groups; ++g)
    {
        LibrpaOrbitalGroup group;
        librpa_orbital_group(result, g, &group);
        if (group.type_index < 0 || group.type_index >= static_cast<int>(species.size()))
            throw std::runtime_error("orbital group species is outside crpa_species_labels");
        groups << g << ' ' << group.atom_index << ' ' << species[group.type_index] << ' ' << group.l
               << ' ' << group.channel << '\n';
    }
    auto weights = output_file(directory, "crpa_orbital_weights.dat");
    weights << "# spin k band group energy_absolute_ev energy_relative_ev occupation kweight "
               "S_norm Mulliken Lowdin\n";
    for (int spin = 0; spin < info.n_spins; ++spin)
        for (int k = 0; k < info.n_kpoints; ++k)
            for (int band = 0; band < info.n_bands; ++band)
            {
                LibrpaOrbitalState state;
                librpa_orbital_state(result, spin, k, band, &state);
                for (int g = 0; g < info.n_groups; ++g)
                    weights << spin << ' ' << k << ' ' << band << ' ' << g << ' '
                            << state.energy_ha * HA2EV << ' '
                            << (state.energy_ha - reference) * HA2EV << ' ' << state.occupation
                            << ' ' << state.kweight << ' ' << state.norm << ' '
                            << state.mulliken_group[g] << ' ' << state.lowdin_group[g] << '\n';
            }
    const int count =
        static_cast<int>(std::floor((params.crpa_pdos_max_ev - params.crpa_pdos_min_ev) /
                                        params.crpa_pdos_step_ev +
                                    1e-10)) +
        1;
    std::vector<double> grid(count);
    for (int e = 0; e < count; ++e)
        grid[e] = reference + (params.crpa_pdos_min_ev + e * params.crpa_pdos_step_ev) / HA2EV;
    const auto length = static_cast<std::size_t>(info.n_spins) * count * info.n_groups;
    std::vector<double> mulliken(length), lowdin(length);
    librpa_orbital_pdos(result, count, grid.data(), params.crpa_pdos_sigma_ev / HA2EV,
                        mulliken.data(), lowdin.data());
    auto pdos = output_file(directory, "crpa_pdos.dat");
    pdos << "# spin energy_relative_ev group Mulliken Lowdin [states/(eV cell "
            "stored-spin-channel)]\n";
    for (int spin = 0; spin < info.n_spins; ++spin)
        for (int e = 0; e < count; ++e)
            for (int g = 0; g < info.n_groups; ++g)
            {
                const auto offset =
                    (static_cast<std::size_t>(spin) * count + e) * info.n_groups + g;
                pdos << spin << ' ' << params.crpa_pdos_min_ev + e * params.crpa_pdos_step_ev << ' '
                     << g << ' ' << mulliken[offset] / HA2EV << ' ' << lowdin[offset] / HA2EV
                     << '\n';
            }
}

void write_result(const LibrpaCrpaResult* result, const std::string& directory, int state_index,
                  int gw_update, bool final_state)
{
    LibrpaCrpaModelInfo info;
    librpa_crpa_model_info(result, &info);
    LibrpaCrpaTensor first;
    librpa_crpa_result_tensor(result, 0, &first);
    const bool has_full_w = first.fully_screened_ri != nullptr;
    auto metadata = output_file(directory, driver::driver_params.crpa_output_format == "dat"
                                               ? "crpa_metadata.dat"
                                               : "crpa_metadata.txt");
    if (state_index >= 0)
        metadata << "state_index = " << state_index << '\n'
                 << "gw_update = " << gw_update << '\n'
                 << "final_state = " << (final_state ? "true" : "false") << '\n'
                 << "response_source = "
                 << (final_state &&
                             driver::get_task(driver::driver_params.task) == driver::task_t::QSGW
                         ? "final_state_response_only"
                         : "shared_entering_state_GW_P0")
                 << '\n';
    metadata << "full_W_computed = " << (has_full_w ? "true" : "false") << '\n';
    metadata
        << "Current loaded eigenstate poles; P_d weights both endpoints by the emitted "
           "probabilities.\n"
        << (info.n_spinor == 2
                ? "Output Phi uses independent B and one jointly normalized spinor frame.\n"
                  "Spinor trial labels are not conserved spin quantum numbers.\n"
                  "Full-spinor J includes all distinct columns; it is not a scalar-shell Hund "
                  "mean.\n"
                : "Output Phi uses independent B window and separate physical spin frames.\n")
        << "First radial shell; eg/t2g labels refer to global Cartesian axes.\n"
        << "Native positive minimax nodes only; no direct U(0) was calculated.\n"
        << "cRPA bare and dielectric/Wr screening use full Coulomb independently of GW flags.\n"
        << (driver::opts.replace_w_head == LIBRPA_SWITCH_ON
                ? (has_full_w ? "Analytic full/constrained head and wings; normalized angular "
                                "Gamma-cell weights.\n"
                              : "Analytic constrained head and wings for U; normalized angular "
                                "Gamma-cell weights.\n")
                : "No optical correction: finite-mesh periodic diagnostic.\n")
        << "Interaction and frequency columns use eV; indices are zero based.\n"
        << "U_intra=mean U_mm,mm; U_inter=mean_(m!=n) U_mm,nn; J=mean U_mn,nm.\n"
        << driver::driver_params.format();
    auto windows = output_file(directory, "crpa_windows.dat");
    auto response_weights = output_file(directory, "crpa_response_weights.dat");
    response_weights << "# spin k band endpoint_weight; Pd multiplies both endpoint weights\n";
    for (int spin = 0; spin < info.n_spins; ++spin)
        for (int k = 0; k < info.n_kpoints; ++k)
        {
            int count;
            const double* weights;
            librpa_crpa_result_weights(result, spin, k, &count, &weights);
            for (int band = 0; band < count; ++band)
                response_weights << spin << ' ' << k << ' ' << band << ' ' << weights[band] << '\n';
        }
    if (driver::driver_params.crpa_response_mode == "orbital_weights")
        windows << "# A is fractional: see crpa_response_weights.dat; no binary A band list\n";
    windows << (state_index >= 0 ? "# window spin k count current_QSGW_band_indices_zero_based\n"
                                 : "# window spin k count original_KS_band_indices_zero_based\n");
    for (int response : {1, 0})
        for (int spin = 0; spin < info.n_spins; ++spin)
            for (int k = 0; k < info.n_kpoints; ++k)
            {
                int count = 0;
                const int* bands = nullptr;
                librpa_crpa_result_bands(result, response, spin, k, &count, &bands);
                windows << (response ? 'A' : 'B') << ' ' << spin << ' ' << k << ' ' << count;
                for (int i = 0; i < count; ++i) windows << ' ' << bands[i];
                windows << '\n';
            }
    auto tensors = output_file(directory, "crpa_tensors.dat");
    auto orbitals = output_file(directory, "crpa_orbitals.dat");
    orbitals << "# atom local_column trial_label\n";
    std::set<int> written_atoms;
    tensors << "# atom spin_left spin_right nu_eV a b c d V_re V_im U_re U_im";
    if (has_full_w) tensors << " W_re W_im";
    tensors << '\n';
    auto summary = output_file(directory, "crpa_summary.dat");
    summary << "# atom spin_left spin_right n_orbitals nu_eV V_intra V_inter J_bare "
               "U_intra U_inter J_U";
    if (has_full_w) summary << " W_intra W_inter J_W";
    summary << '\n';
    for (int record = 0; record < librpa_crpa_result_size(result); ++record)
    {
        LibrpaCrpaTensor tensor;
        librpa_crpa_result_tensor(result, record, &tensor);
        const int n = tensor.n_orbitals;
        if (written_atoms.insert(tensor.atom_index).second)
            for (int a = 0; a < n; ++a)
                orbitals << tensor.atom_index << ' ' << a << ' '
                         << librpa_crpa_result_orbital_label(result, tensor.atom_index, a) << '\n';
        const double frequency = tensor.frequency_ha * librpa_int::HA2EV;
        summary << tensor.atom_index << ' ' << tensor.spin_left << ' ' << tensor.spin_right << ' '
                << n << ' ' << frequency;
        for (const auto* values :
             {tensor.bare_ri, tensor.partially_screened_ri, tensor.fully_screened_ri})
            if (values)
            {
                std::array<double, 6> averages;
                librpa_crpa_tensor_averages(n, values, averages.data());
                for (int channel = 0; channel < 3; ++channel)
                    summary << ' ' << averages[2 * channel] * librpa_int::HA2EV;
            }
        summary << '\n';
        for (int a = 0; a < n; ++a)
            for (int b = 0; b < n; ++b)
                for (int c = 0; c < n; ++c)
                    for (int d = 0; d < n; ++d)
                    {
                        const int index = 2 * ((a * n + b) * n * n + c * n + d);
                        tensors << tensor.atom_index << ' ' << tensor.spin_left << ' '
                                << tensor.spin_right << ' ' << frequency << ' ' << a << ' ' << b
                                << ' ' << c << ' ' << d;
                        for (const auto* values : {tensor.bare_ri, tensor.partially_screened_ri,
                                                   tensor.fully_screened_ri})
                            if (values)
                                tensors << ' ' << values[index] * librpa_int::HA2EV << ' '
                                        << values[index + 1] * librpa_int::HA2EV;
                        tensors << '\n';
                    }
    }
}

}  // namespace

namespace driver
{
namespace
{
std::string crpa_tracking_definition()
{
    const auto& p = driver_params;
    std::ostringstream text;
    text << std::setprecision(17) << p.crpa_species_labels << '|' << p.crpa_correlated_species
         << '|' << p.crpa_ligand_species << '|' << p.crpa_parent_orbitals << '|'
         << p.crpa_output_orbitals << '|' << p.crpa_response_mode << '|'
         << p.crpa_tracking_min_capture;
    for (const double edge : p.crpa_orbital_windows_ha) text << '|' << edge;
    text << ";bands";
    for (const int band : p.crpa_orbital_bands) text << '|' << band;
    return text.str();
}
}  // namespace

std::string CrpaDriverInput::tracking_definition() const { return crpa_tracking_definition(); }

CrpaDriverInput read_crpa_driver_input()
{
    using librpa_int::global::mpi_comm_global_h;
    CrpaDriverInput result;
    auto& species = result.species;
    auto& overlap_packed = result.overlap_packed;
    std::string error;
    try
    {
        int n_spins, n_bands;
        librpa_get_meanfield_dimensions(h.get_c_handler(), 0, &n_spins, &result.n_kpoints, &n_bands,
                                        &result.n_aos, &result.n_spinor);
        std::string species_text = driver_params.crpa_species_labels;
        std::replace(species_text.begin(), species_text.end(), ',', ' ');
        std::istringstream stream(species_text);
        std::string label;
        while (stream >> label) species.push_back(label);
        const auto filename =
            librpa_int::is_absolute_path(driver_params.crpa_overlap_file)
                ? driver_params.crpa_overlap_file
                : librpa_int::join_path(driver_params.input_dir, driver_params.crpa_overlap_file);
        const auto overlap = librpa::reader::read_abacus_nao_overlap_csr(filename);
        std::vector<double> kfrac(3 * result.n_kpoints);
        librpa_get_crpa_kgrid(h.get_c_handler(), result.n_kpoints, kfrac.data());
        std::vector<librpa_int::Vector3_Order<double>> kpoints;
        for (int k = 0; k < result.n_kpoints; ++k)
            kpoints.push_back({kfrac[3 * k], kfrac[3 * k + 1], kfrac[3 * k + 2]});
        const auto overlap_k = librpa::reader::fourier_nao_overlap(overlap, kpoints);
        for (const auto& matrix : overlap_k)
        {
            const auto spatial =
                librpa::reader::spatial_nao_overlap(matrix, result.n_aos, result.n_spinor);
            for (int i = 0; i < spatial.size; ++i)
            {
                overlap_packed.push_back(spatial.c[i].real());
                overlap_packed.push_back(spatial.c[i].imag());
            }
        }
    }
    catch (const std::exception& failure)
    {
        error = failure.what();
    }
    const int failed = !error.empty();
    int any_failed = 0;
    mpi_comm_global_h.allreduce(&failed, &any_failed, 1, MPI_MAX);
    if (any_failed)
        throw std::runtime_error("cRPA native input: " +
                                 (error.empty() ? "another MPI rank failed" : error));

    return result;
}

LibrpaCrpaInput CrpaDriverInput::input()
{
    LibrpaCrpaInput input;
    librpa_init_crpa_input(&input);
    labels.clear();
    for (const auto& name : species) labels.push_back(name.c_str());
    input.n_species = static_cast<int>(labels.size());
    input.species_labels = labels.data();
    input.correlated_species = driver_params.crpa_correlated_species.c_str();
    input.ligand_species = driver_params.crpa_ligand_species.c_str();
    input.parent_orbitals = driver_params.crpa_parent_orbitals.c_str();
    input.output_orbitals = driver_params.crpa_output_orbitals.c_str();
    input.use_orbital_weights = driver_params.crpa_response_mode == "orbital_weights"
                                    ? LIBRPA_SWITCH_ON
                                    : LIBRPA_SWITCH_OFF;
    input.n_response_edges = static_cast<int>(driver_params.crpa_response_windows_ha.size());
    input.response_windows_ha = driver_params.crpa_response_windows_ha.data();
    input.n_orbital_edges = static_cast<int>(driver_params.crpa_orbital_windows_ha.size());
    input.orbital_windows_ha = driver_params.crpa_orbital_windows_ha.data();
    input.n_response_bands = static_cast<int>(driver_params.crpa_response_bands.size());
    input.response_bands = driver_params.crpa_response_bands.data();
    input.n_orbital_bands = static_cast<int>(driver_params.crpa_orbital_bands.size());
    input.orbital_bands = driver_params.crpa_orbital_bands.data();
    input.n_kpoints = n_kpoints;
    input.n_aos = n_aos;
    input.overlap_k_ri = overlap_packed.data();
    input.residual_tol = driver_params.crpa_residual_tol;
    input.retain_raw_model = driver_params.crpa_output_raw ? LIBRPA_SWITCH_ON : LIBRPA_SWITCH_OFF;
    input.tracking_min_capture = driver_params.crpa_tracking_min_capture;
    return input;
}

void write_crpa_state_outputs(const LibrpaCrpaResult* result, CrpaDriverInput& input,
                              const std::string& directory, int state_index, int gw_update,
                              bool final_state, bool write_analysis,
                              const LibrpaOrbitalResult* supplied_analysis)
{
    const auto& comm = librpa_int::global::mpi_comm_global_h;
    std::string error;
    try
    {
        librpa_int::create_directories(directory.c_str(), comm.myid);
    }
    catch (const std::exception& failure)
    {
        error = failure.what();
    }
    int failed = !error.empty(), any_failed = 0;
    comm.allreduce(&failed, &any_failed, 1, MPI_MAX);
    if (any_failed) throw std::runtime_error("cRPA output directory: " + error);
    comm.barrier();

    librpa::OrbitalResult owned_analysis(nullptr, librpa_delete_orbital_result);
    if (!supplied_analysis && write_analysis && driver_params.crpa_output_analysis)
        owned_analysis = h.analyze_orbitals(input.n_kpoints, input.n_aos, input.overlap_packed,
                                            driver_params.crpa_residual_tol);
    const auto* analysis = supplied_analysis ? supplied_analysis : owned_analysis.get();
    error.clear();
    try
    {
        if (comm.is_root())
        {
            if (analysis) write_orbital_analysis(analysis, input.species, directory);
            if (result)
            {
                write_result(result, directory, state_index, gw_update, final_state);
                write_crpa_model(result, directory, state_index, gw_update, final_state);
                if (driver_params.crpa_track_orbitals)
                    librpa_write_crpa_tracking_checkpoint(
                        result, librpa_int::join_path(directory, "crpa_tracking.dat").c_str(),
                        state_index, input.tracking_definition().c_str(), 0, nullptr);
            }
        }
    }
    catch (const std::exception& failure)
    {
        error = failure.what();
    }
    failed = !error.empty();
    comm.allreduce(&failed, &any_failed, 1, MPI_MAX);
    if (any_failed)
        throw std::runtime_error("cRPA output: " +
                                 (error.empty() ? "another MPI rank failed" : error));
}

void task_crpa()
{
    auto config = read_crpa_driver_input();
    write_crpa_state_outputs(nullptr, config, opts.output_dir);
    if (driver_params.crpa_analysis_only)
    {
        return;
    }
    if (opts.replace_w_head == LIBRPA_SWITCH_ON) prepare_crpa_optical_velocity();
    auto input = config.input();
    librpa::CrpaResult result(librpa_compute_crpa_window(h.get_c_handler(), &opts, &input),
                              librpa_delete_crpa_result);
    write_crpa_state_outputs(result.get(), config, opts.output_dir, -1, -1, false, false);
    librpa_int::global::lib_printf_root("cRPA tensors and orbital analysis written to %s\n",
                                        opts.output_dir);
}
}  // namespace driver
