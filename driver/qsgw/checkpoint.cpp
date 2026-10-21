#include "checkpoint.h"

#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <stdexcept>

namespace librpa_int
{
namespace qsgw
{
namespace
{
constexpr char CHECKPOINT_MAGIC[] = "LIBRPA_QSGW";
constexpr std::uint32_t CHECKPOINT_BYTE_ORDER = 0x01020304;
// ASCII format tags, not content hashes. Preserve the existing checkpoint bytes.
constexpr std::uint32_t REFERENCE_VELOCITY_TAG = 0x48575631;  // HWV1
constexpr std::uint32_t NATIVE_VELOCITY_TAG = 0x48575632;     // HWV2

template <typename T>
void put(std::ostream& output, const T& value)
{
    output.write(reinterpret_cast<const char*>(&value), sizeof(T));
}

template <typename T>
T get(std::istream& input)
{
    T value{};
    input.read(reinterpret_cast<char*>(&value), sizeof(T));
    if (!input) throw std::runtime_error("Truncated QSGW checkpoint");
    return value;
}

double get_real(std::istream& input)
{
    const double value = get<double>(input);
    if (!std::isfinite(value)) throw std::runtime_error("Non-finite QSGW checkpoint value");
    return value;
}

void put_complex(std::ostream& output, const cplxdb value)
{
    put(output, value.real());
    put(output, value.imag());
}

cplxdb get_complex(std::istream& input)
{
    const double real = get_real(input);
    const double imaginary = get_real(input);
    return {real, imaginary};
}

void require_equal(const double saved, const double expected)
{
    if (!std::isfinite(expected) || std::abs(saved - expected) > 1.0e-12)
        throw std::runtime_error(
            "QSGW checkpoint reference or update settings differ from the input");
}

void write_meanfield(std::ostream& output, const MeanField& meanfield)
{
    for (const int count :
         {meanfield.get_n_spins(), meanfield.get_n_kpoints(), meanfield.get_n_bands(),
          meanfield.get_n_aos(), meanfield.get_n_spinor()})
        put(output, static_cast<std::int32_t>(count));
    put(output, meanfield.get_efermi());
    for (int spin = 0; spin < meanfield.get_n_spins(); ++spin)
    {
        for (int k = 0; k < meanfield.get_n_kpoints(); ++k)
            for (int band = 0; band < meanfield.get_n_bands(); ++band)
            {
                put(output, meanfield.get_eigenvals()[spin](k, band));
                put(output, meanfield.get_weight()[spin](k, band));
            }
        for (int spinor = 0; spinor < meanfield.get_n_spinor(); ++spinor)
            for (int k = 0; k < meanfield.get_n_kpoints(); ++k)
            {
                const auto& wfc = meanfield.get_eigenvectors().at(spin).at(spinor).at(k);
                if (wfc.nr != meanfield.get_n_bands() || wfc.nc != meanfield.get_n_aos())
                    throw std::runtime_error("QSGW checkpoint requires replicated wavefunctions");
                for (int i = 0; i < wfc.size; ++i) put_complex(output, wfc.c[i]);
            }
    }
}

// Shapes are checked against native inputs before any checkpoint allocation.
MeanField read_meanfield(std::istream& input, const MeanField& reference)
{
    for (const int count :
         {reference.get_n_spins(), reference.get_n_kpoints(), reference.get_n_bands(),
          reference.get_n_aos(), reference.get_n_spinor()})
        if (get<std::int32_t>(input) != count)
            throw std::runtime_error("QSGW checkpoint mean-field dimensions differ from the input");
    MeanField result = reference;
    result.get_efermi() = get_real(input);
    for (int spin = 0; spin < result.get_n_spins(); ++spin)
    {
        for (int k = 0; k < result.get_n_kpoints(); ++k)
            for (int band = 0; band < result.get_n_bands(); ++band)
            {
                result.get_eigenvals()[spin](k, band) = get_real(input);
                result.get_weight()[spin](k, band) = get_real(input);
            }
        for (int spinor = 0; spinor < result.get_n_spinor(); ++spinor)
            for (int k = 0; k < result.get_n_kpoints(); ++k)
            {
                auto& wfc = result.get_eigenvectors().at(spin).at(spinor).at(k);
                for (int i = 0; i < wfc.size; ++i) wfc.c[i] = get_complex(input);
            }
    }
    return result;
}

void compare_meanfield(const MeanField& saved, const MeanField& reference)
{
    require_equal(saved.get_efermi(), reference.get_efermi());
    for (int spin = 0; spin < reference.get_n_spins(); ++spin)
    {
        for (int k = 0; k < reference.get_n_kpoints(); ++k)
            for (int band = 0; band < reference.get_n_bands(); ++band)
            {
                require_equal(saved.get_eigenvals()[spin](k, band),
                              reference.get_eigenvals()[spin](k, band));
                require_equal(saved.get_weight()[spin](k, band),
                              reference.get_weight()[spin](k, band));
            }
        for (int spinor = 0; spinor < reference.get_n_spinor(); ++spinor)
            for (int k = 0; k < reference.get_n_kpoints(); ++k)
            {
                const auto& expected = reference.get_eigenvectors().at(spin).at(spinor).at(k);
                const auto& actual = saved.get_eigenvectors().at(spin).at(spinor).at(k);
                for (int i = 0; i < expected.size; ++i)
                {
                    require_equal(actual.c[i].real(), expected.c[i].real());
                    require_equal(actual.c[i].imag(), expected.c[i].imag());
                }
            }
    }
}

void write_matrices(std::ostream& output, const SpinKMatrixMap& matrices,
                    const MeanField& reference)
{
    const int n = reference.get_n_bands();
    for (int spin = 0; spin < reference.get_n_spins(); ++spin)
        for (int k = 0; k < reference.get_n_kpoints(); ++k)
        {
            const auto& matrix = matrices.at(spin).at(k);
            if (matrix.nr() != n || matrix.nc() != n)
                throw std::runtime_error("Invalid QSGW checkpoint Hamiltonian dimensions");
            for (int row = 0; row < n; ++row)
                for (int col = 0; col < n; ++col) put_complex(output, matrix(row, col));
        }
}

SpinKMatrixMap read_matrices(std::istream& input, const MeanField& reference)
{
    SpinKMatrixMap result;
    const int n = reference.get_n_bands();
    for (int spin = 0; spin < reference.get_n_spins(); ++spin)
        for (int k = 0; k < reference.get_n_kpoints(); ++k)
        {
            Matz matrix(n, n, MAJOR::ROW);
            for (int row = 0; row < n; ++row)
                for (int col = 0; col < n; ++col) matrix(row, col) = get_complex(input);
            result[spin].emplace(k, std::move(matrix));
        }
    return result;
}

bool live_velocity_checkpoint(const CheckpointContext& context)
{
    return context.head_mode != 0 && context.velocity_mode != "reference";
}

bool native_velocity_checkpoint(const CheckpointContext& context)
{
    return live_velocity_checkpoint(context) && context.velocity_mode == "current_hr";
}

bool shifted_cut_checkpoint(const CheckpointContext& context)
{
    if (context.cut_mode != 0 && context.cut_mode != 1)
        throw std::runtime_error("QSGW checkpoint cut mode must be 0 or 1");
    return context.cut_mode == 1;
}

int velocity_mode_code(const CheckpointContext& context)
{
    if (context.velocity_mode == "current_hr") return 1;
    if (context.velocity_mode == "rotate_reference") return 2;
    throw std::runtime_error("Unsupported QSGW checkpoint velocity convention");
}

void require_native_velocity_context(const CheckpointContext& context)
{
    if (!context.ao_operators || !context.velocity_period || !context.velocity_lattice)
        throw std::runtime_error("Missing current-HR checkpoint operators or geometry");
    const auto& operators = *context.ao_operators;
    if (operators.dimension != context.reference.get_n_aos() * context.reference.get_n_spinor() ||
        operators.native_dft_hamiltonian.size() !=
            static_cast<std::size_t>(context.reference.get_n_spins()))
        throw std::runtime_error(
            "QSGW current-HR checkpoint requires one native DFT Hamiltonian per spin channel");
    for (const auto& blocks : operators.native_dft_hamiltonian)
    {
        if (blocks.empty())
            throw std::runtime_error(
                "Missing native DFT Hamiltonian in QSGW current-HR checkpoint");
        for (const auto& [r, matrix] : blocks)
        {
            if (matrix.nr() != operators.dimension || matrix.nc() != operators.dimension)
                throw std::runtime_error("Invalid native DFT checkpoint Hamiltonian dimensions");
            for (int i = 0; i < matrix.nr(); ++i)
                for (int j = 0; j < matrix.nc(); ++j)
                    if (!std::isfinite(matrix(i, j).real()) || !std::isfinite(matrix(i, j).imag()))
                        throw std::runtime_error("Non-finite native DFT checkpoint Hamiltonian");
        }
    }
}

// Save the immutable numerical operators, not their filenames. A resumed run
// must use the same map even if files move or keep their names but change.
void write_velocity_context(std::ostream& output, const CheckpointContext& context)
{
    const int mode = velocity_mode_code(context);
    // HWV2 binds current-HR velocities to the method-2 native DFT baseline.
    // Reference rotation retains its existing HWV1 numerical contract.
    put(output, mode == 1 ? NATIVE_VELOCITY_TAG : REFERENCE_VELOCITY_TAG);
    put(output, std::int32_t(mode));
    if (mode == 1)
    {
        require_native_velocity_context(context);
        const auto& period = *context.velocity_period;
        for (int value : {period.x, period.y, period.z}) put(output, std::int32_t(value));
        const auto& lattice = *context.velocity_lattice;
        for (double value : {lattice.e11, lattice.e12, lattice.e13, lattice.e21, lattice.e22,
                             lattice.e23, lattice.e31, lattice.e32, lattice.e33})
            put(output, value);
        put(output, std::int32_t(context.ao_operators->dimension));
        const auto write_localized = [&](const LocalizedMatrix& blocks)
        {
            put(output, std::uint64_t(blocks.size()));
            for (const auto& [r, matrix] : blocks)
            {
                for (int value : {r.x, r.y, r.z}) put(output, std::int32_t(value));
                for (int i = 0; i < matrix.nr(); ++i)
                    for (int j = 0; j < matrix.nc(); ++j) put_complex(output, matrix(i, j));
            }
        };
        write_localized(context.ao_operators->overlap);
        for (const auto& blocks : context.ao_operators->position) write_localized(blocks);
        const auto& native_h = context.ao_operators->native_dft_hamiltonian;
        put(output, std::uint32_t(native_h.size()));
        for (const auto& blocks : native_h) write_localized(blocks);
    }
    else
    {
        if (!context.reference_velocity)
            throw std::runtime_error("Missing reference velocity for QSGW checkpoint");
        for (const auto& by_k : *context.reference_velocity)
            for (const auto& by_direction : by_k)
                for (const auto& matrix : by_direction)
                    for (int i = 0; i < matrix.size; ++i) put_complex(output, matrix.c[i]);
    }
}

void check_velocity_context(std::istream& input, const CheckpointContext& context)
{
    const int mode = velocity_mode_code(context);
    const auto expected_marker = mode == 1 ? NATIVE_VELOCITY_TAG : REFERENCE_VELOCITY_TAG;
    if (get<std::uint32_t>(input) != expected_marker || get<std::int32_t>(input) != mode)
        throw std::runtime_error(
            "QSGW checkpoint head/wing velocity convention differs from the input");
    if (context.velocity_mode == "current_hr")
    {
        require_native_velocity_context(context);
        const auto& period = *context.velocity_period;
        for (int value : {period.x, period.y, period.z})
            require_equal(get<std::int32_t>(input), value);
        const auto& lattice = *context.velocity_lattice;
        for (double value : {lattice.e11, lattice.e12, lattice.e13, lattice.e21, lattice.e22,
                             lattice.e23, lattice.e31, lattice.e32, lattice.e33})
            require_equal(get_real(input), value);
        const int n = context.ao_operators->dimension;
        require_equal(get<std::int32_t>(input), n);
        const auto check_localized = [&](const LocalizedMatrix& blocks)
        {
            if (get<std::uint64_t>(input) != blocks.size())
                throw std::runtime_error("QSGW checkpoint AO support differs from the input");
            for (const auto& [r, matrix] : blocks)
            {
                for (int value : {r.x, r.y, r.z}) require_equal(get<std::int32_t>(input), value);
                for (int i = 0; i < n; ++i)
                    for (int j = 0; j < n; ++j)
                    {
                        const auto value = get_complex(input);
                        require_equal(value.real(), matrix(i, j).real());
                        require_equal(value.imag(), matrix(i, j).imag());
                    }
            }
        };
        check_localized(context.ao_operators->overlap);
        for (const auto& blocks : context.ao_operators->position) check_localized(blocks);
        const auto& native_h = context.ao_operators->native_dft_hamiltonian;
        if (get<std::uint32_t>(input) != native_h.size())
            throw std::runtime_error(
                "QSGW checkpoint native DFT spin channels differ from the input");
        for (const auto& blocks : native_h) check_localized(blocks);
    }
    else
    {
        if (!context.reference_velocity)
            throw std::runtime_error("Missing reference velocity for QSGW checkpoint");
        for (const auto& by_k : *context.reference_velocity)
            for (const auto& by_direction : by_k)
                for (const auto& matrix : by_direction)
                    for (int i = 0; i < matrix.size; ++i)
                    {
                        const auto value = get_complex(input);
                        require_equal(value.real(), matrix.c[i].real());
                        require_equal(value.imag(), matrix.c[i].imag());
                    }
    }
}

void write_context(std::ostream& output, const CheckpointContext& context)
{
    put(output, std::int32_t(context.band_reference != nullptr));
    for (const int value :
         {context.cut_mode, context.unoccupied_keep, context.nfreq, context.head_mode})
        put(output, static_cast<std::int32_t>(value));
    put(output, context.cut_shift_ha);
    put(output, context.mixing_beta);
    put(output, std::int32_t(0));  // Existing fixed-mixer code in the binary format.
    put(output, std::int32_t(0));  // Existing unweighted-residual metric code.
    for (const auto* points : {&context.kpoints, &context.band_kpoints})
    {
        put(output, static_cast<std::int32_t>(points->size()));
        for (const auto& point : *points)
            for (const double value : {point.x, point.y, point.z}) put(output, value);
    }
    put(output, static_cast<std::int32_t>(context.k_weights.size()));
    for (const double weight : context.k_weights) put(output, weight);
    write_meanfield(output, context.reference);
    write_matrices(output, context.dft_vxc, context.reference);
    if (context.band_reference)
    {
        write_meanfield(output, *context.band_reference);
        write_matrices(output, context.dft_vxc_band, *context.band_reference);
    }
}

void check_context(std::istream& input, const CheckpointContext& context,
                   const std::uint32_t version)
{
    if (get<std::int32_t>(input) != (context.band_reference != nullptr))
        throw std::runtime_error("QSGW checkpoint task differs from the input");
    for (const int value :
         {context.cut_mode, context.unoccupied_keep, context.nfreq, context.head_mode})
        require_equal(get<std::int32_t>(input), value);
    require_equal(get_real(input), context.cut_shift_ha);
    require_equal(get_real(input), context.mixing_beta);
    // Fixed-mode files retain compatibility with version 1--3 checkpoints.
    // Mixer history belongs to the separate automatic-mixing implementation.
    if (version >= 2 && get<std::int32_t>(input) != 0)
        throw std::runtime_error("Automatic-mixing checkpoints require a mixer-aware build");
    if (version >= 3 && get<std::int32_t>(input) != 0)
        throw std::runtime_error("QSGW checkpoint requires an unsupported mixing metric");
    for (const auto* points : {&context.kpoints, &context.band_kpoints})
    {
        if (get<std::int32_t>(input) != static_cast<int>(points->size()))
            throw std::runtime_error("QSGW checkpoint k-point count differs from the input");
        for (const auto& point : *points)
            for (const double value : {point.x, point.y, point.z})
                require_equal(get_real(input), value);
    }
    if (get<std::int32_t>(input) != static_cast<int>(context.k_weights.size()))
        throw std::runtime_error("QSGW checkpoint k-weight count differs from the input");
    for (const double weight : context.k_weights) require_equal(get_real(input), weight);
    const auto check_channel = [&](const MeanField& reference, const SpinKMatrixMap& vxc)
    {
        compare_meanfield(read_meanfield(input, reference), reference);
        const auto saved_vxc = read_matrices(input, reference);
        for (const auto& [spin, by_k] : saved_vxc)
            for (const auto& [k, matrix] : by_k)
                for (int row = 0; row < matrix.nr(); ++row)
                    for (int col = 0; col < matrix.nc(); ++col)
                    {
                        const auto expected = vxc.at(spin).at(k)(row, col);
                        require_equal(matrix(row, col).real(), expected.real());
                        require_equal(matrix(row, col).imag(), expected.imag());
                    }
    };
    check_channel(context.reference, context.dft_vxc);
    if (context.band_reference) check_channel(*context.band_reference, context.dft_vxc_band);
}

}  // namespace

void write_qsgw_checkpoint(const std::string& filename, const CheckpointContext& context,
                           const IterationSummary& summary, const MeanField& live,
                           const MeanField* band_live, const SpinKMatrixMap& hamiltonian,
                           const SpinKMatrixMap& band_hamiltonian)
{
    const std::string temporary = filename + ".tmp";
    std::ofstream output(temporary, std::ios::binary | std::ios::trunc);
    if (!output) throw std::runtime_error("Cannot write QSGW checkpoint: " + temporary);
    output.write(CHECKPOINT_MAGIC, sizeof(CHECKPOINT_MAGIC) - 1);
    // Versions 8/9 retain the v3/v5 payload and identify shifted cut mode 1
    // on either grid-only or band-enabled runs. Earlier mode 1 was unshifted;
    // versions 6/7 identified the older grid-only mode-2 shifted cut.
    // Versions 10/11 replace whole-H current-HR reconstruction with native
    // DFT plus correction interpolation and bind H0(R) in the velocity context.
    // Version 4 belongs to another equation/mixer branch. Do not reinterpret it.
    const bool live_velocity = live_velocity_checkpoint(context);
    const bool native_velocity = native_velocity_checkpoint(context);
    const bool shifted_cut = shifted_cut_checkpoint(context);
    const std::uint32_t version =
        native_velocity ? (shifted_cut ? 11 : 10)
                        : (shifted_cut ? (live_velocity ? 9 : 8) : (live_velocity ? 5 : 3));
    put(output, version);
    put(output, CHECKPOINT_BYTE_ORDER);
    write_context(output, context);
    if (live_velocity) write_velocity_context(output, context);
    put(output, std::int32_t(summary.iteration));
    for (const double value :
         {summary.maximum_eigenvalue_change_ev, summary.residual_l2_ha, summary.residual_max_ha,
          summary.fermi_energy_ev, summary.gap_ev, summary.electron_count, summary.beta})
        put(output, value);
    put(output, std::int32_t(summary.converged));
    write_meanfield(output, live);
    write_matrices(output, hamiltonian, context.reference);
    if (context.band_reference)
    {
        write_meanfield(output, *band_live);
        write_matrices(output, band_hamiltonian, *context.band_reference);
    }
    output.flush();
    if (!output) throw std::runtime_error("Failed to flush QSGW checkpoint: " + temporary);
    output.close();
    if (!output) throw std::runtime_error("Failed to close QSGW checkpoint: " + temporary);
    // Keep the last complete iteration if the process dies during a later write.
    std::filesystem::rename(temporary, filename);
}

IterationSummary read_qsgw_checkpoint(const std::string& filename, const CheckpointContext& context,
                                      MeanField& live, MeanField* band_live,
                                      SpinKMatrixMap& hamiltonian, SpinKMatrixMap& band_hamiltonian)
{
    std::ifstream input(filename, std::ios::binary);
    if (!input) throw std::runtime_error("Cannot read QSGW checkpoint: " + filename);
    char magic[sizeof(CHECKPOINT_MAGIC) - 1];
    input.read(magic, sizeof(magic));
    if (!input || std::string(magic, sizeof(magic)) != CHECKPOINT_MAGIC)
        throw std::runtime_error("Unsupported QSGW checkpoint format");
    const auto version = get<std::uint32_t>(input);
    if ((version < 1 || (version > 3 && version != 5 && version != 6 && version != 7 &&
                         version != 8 && version != 9 && version != 10 && version != 11)) ||
        get<std::uint32_t>(input) != CHECKPOINT_BYTE_ORDER)
        throw std::runtime_error("Unsupported QSGW checkpoint format or byte order");
    if (version == 6 || version == 7)
        throw std::runtime_error(
            "QSGW checkpoint versions 6/7 use the historical mode-2 cut; "
            "automatic migration to shifted mode 1 is unsupported");
    const bool native_velocity = version == 10 || version == 11;
    const bool requested_native_velocity = native_velocity_checkpoint(context);
    if (requested_native_velocity && !native_velocity)
        throw std::runtime_error(
            "QSGW current-HR native-baseline velocities require checkpoint version 10 or 11; "
            "legacy whole-H velocity checkpoints cannot continue with this convention");
    if (native_velocity && !requested_native_velocity)
        throw std::runtime_error(
            "QSGW checkpoint uses a different native-baseline head/wing velocity convention");
    const bool shifted_cut = version == 8 || version == 9 || version == 11;
    const bool requested_shifted_cut = shifted_cut_checkpoint(context);
    if (requested_shifted_cut && !shifted_cut)
        throw std::runtime_error(
            "Legacy QSGW cut checkpoints cannot continue shifted mode 1; "
            "this convention requires checkpoint version 8, 9 or 11");
    if (shifted_cut != requested_shifted_cut)
        throw std::runtime_error("QSGW checkpoint shifted-cut convention differs from the input");
    const bool live_velocity = version == 5 || version == 9 || native_velocity;
    if (live_velocity != live_velocity_checkpoint(context))
        throw std::runtime_error("QSGW checkpoint lacks the requested head/wing velocity contract");
    check_context(input, context, version);
    if (live_velocity) check_velocity_context(input, context);
    IterationSummary summary;
    summary.iteration = get<std::int32_t>(input);
    if (summary.iteration < 1) throw std::runtime_error("Invalid QSGW checkpoint iteration");
    summary.maximum_eigenvalue_change_ev = get_real(input);
    summary.residual_l2_ha = get_real(input);
    summary.residual_max_ha = get_real(input);
    summary.fermi_energy_ev = get_real(input);
    summary.gap_ev = get_real(input);
    summary.electron_count = get_real(input);
    summary.beta = get_real(input);
    const int converged = get<std::int32_t>(input);
    if (converged != 0 && converged != 1)
        throw std::runtime_error("Invalid QSGW checkpoint convergence flag");
    summary.converged = converged != 0;
    MeanField restored = read_meanfield(input, context.reference);
    SpinKMatrixMap restored_h = read_matrices(input, context.reference);
    MeanField restored_band;
    SpinKMatrixMap restored_band_h;
    if (context.band_reference)
    {
        restored_band = read_meanfield(input, *context.band_reference);
        restored_band_h = read_matrices(input, *context.band_reference);
    }
    if (input.peek() != std::char_traits<char>::eof())
        throw std::runtime_error("Unexpected trailing QSGW checkpoint data");
    live = std::move(restored);
    hamiltonian = std::move(restored_h);
    if (context.band_reference)
    {
        *band_live = std::move(restored_band);
        band_hamiltonian = std::move(restored_band_h);
    }
    return summary;
}

void broadcast_checkpoint_meanfield(MeanField& live, const MpiCommHandler& communicator)
{
    communicator.bcast(&live.get_efermi(), 1, 0);
    for (int spin = 0; spin < live.get_n_spins(); ++spin)
    {
        auto& eigenvalues = live.get_eigenvals()[spin];
        auto& occupations = live.get_weight()[spin];
        communicator.bcast(eigenvalues.c, eigenvalues.nr * eigenvalues.nc, 0);
        communicator.bcast(occupations.c, occupations.nr * occupations.nc, 0);
        for (int spinor = 0; spinor < live.get_n_spinor(); ++spinor)
            for (int k = 0; k < live.get_n_kpoints(); ++k)
            {
                auto& wfc = live.get_eigenvectors().at(spin).at(spinor).at(k);
                communicator.bcast(wfc.c, wfc.size, 0);
            }
    }
}

}  // namespace qsgw
}  // namespace librpa_int
