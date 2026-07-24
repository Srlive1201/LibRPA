#include "hamiltonian_mixing.h"

#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>

namespace librpa_int
{
namespace qsgw
{
namespace
{

constexpr double hermitian_tolerance = 1.0e-10;

using MatrixBlockLayout = detail::HamiltonianMatrixBlockLayout;
using MapLayout = detail::HamiltonianMapLayout;

MapLayout build_layout(const SpinKMatrixMap& values, const char* channel)
{
    if (values.empty())
    {
        throw std::invalid_argument(
            std::string("QSGW ") + channel + " Hamiltonian map is empty");
    }

    MapLayout layout;
    for (const auto& [spin, by_kpoint] : values)
    {
        if (by_kpoint.empty())
        {
            throw std::invalid_argument(
                std::string("QSGW ") + channel +
                " Hamiltonian contains an empty spin map");
        }
        for (const auto& [kpoint, value] : by_kpoint)
        {
            if (value.nr() <= 0 || value.nr() != value.nc())
            {
                throw std::invalid_argument(
                    std::string("QSGW ") + channel +
                    " Hamiltonian matrix must be non-empty and square");
            }
            for (int row = 0; row < value.nr(); ++row)
            {
                for (int column = 0; column < value.nc(); ++column)
                {
                    const cplxdb element = value(row, column);
                    if (!std::isfinite(element.real()) ||
                        !std::isfinite(element.imag()))
                    {
                        throw std::invalid_argument(
                            std::string("QSGW ") + channel +
                            " Hamiltonian contains non-finite data");
                    }
                    if (std::abs(element - std::conj(value(column, row))) >
                        hermitian_tolerance)
                    {
                        throw std::invalid_argument(
                            std::string("QSGW ") + channel +
                            " Hamiltonian is not Hermitian");
                    }
                }
            }

            MatrixBlockLayout block;
            block.spin = spin;
            block.kpoint = kpoint;
            block.rows = value.nr();
            block.columns = value.nc();
            block.major = value.major();
            block.packed_offset = layout.packed_size;
            layout.packed_size += 2 * value.nr() * value.nc();
            layout.blocks.push_back(block);
        }
    }
    return layout;
}

void require_matching_layout(const SpinKMatrixMap& values,
                             const MapLayout& expected,
                             const char* channel)
{
    const MapLayout actual = build_layout(values, channel);
    if (actual.blocks.size() != expected.blocks.size() ||
        actual.packed_size != expected.packed_size)
    {
        throw std::invalid_argument(
            std::string("QSGW ") + channel +
            " Hamiltonian layout changed during mixing");
    }
    for (std::size_t index = 0; index < expected.blocks.size(); ++index)
    {
        const MatrixBlockLayout& lhs = actual.blocks[index];
        const MatrixBlockLayout& rhs = expected.blocks[index];
        if (lhs.spin != rhs.spin || lhs.kpoint != rhs.kpoint ||
            lhs.rows != rhs.rows || lhs.columns != rhs.columns)
        {
            throw std::invalid_argument(
                std::string("QSGW ") + channel +
                " Hamiltonian keys or dimensions changed during mixing");
        }
    }
}

matrix pack_map(const SpinKMatrixMap& values, const MapLayout& layout)
{
    matrix packed(layout.packed_size, 1, true);
    for (const MatrixBlockLayout& block : layout.blocks)
    {
        const Matz& value = values.at(block.spin).at(block.kpoint);
        int offset = block.packed_offset;
        for (int row = 0; row < block.rows; ++row)
        {
            for (int column = 0; column < block.columns; ++column)
            {
                packed(offset++, 0) = value(row, column).real();
                packed(offset++, 0) = value(row, column).imag();
            }
        }
    }
    return packed;
}

SpinKMatrixMap unpack_map(const matrix& packed, const MapLayout& layout)
{
    if (packed.nr != layout.packed_size || packed.nc != 1)
    {
        throw std::runtime_error(
            "QSGW packed Hamiltonian has an unexpected shape");
    }

    SpinKMatrixMap result;
    for (const MatrixBlockLayout& block : layout.blocks)
    {
        Matz value(block.rows, block.columns, block.major);
        int offset = block.packed_offset;
        for (int row = 0; row < block.rows; ++row)
        {
            for (int column = 0; column < block.columns; ++column)
            {
                value(row, column) =
                    cplxdb(packed(offset, 0), packed(offset + 1, 0));
                offset += 2;
            }
        }
        result[block.spin][block.kpoint] = std::move(value);
    }
    return result;
}

} // namespace

SpinKHamiltonianResidual measure_spin_k_hamiltonian_residual(
    const SpinKMatrixMap& output,
    const SpinKMatrixMap& input)
{
    const MapLayout layout = build_layout(input, "input");
    require_matching_layout(output, layout, "output");
    double square_sum = 0.0;
    double maximum = 0.0;
    for (const MatrixBlockLayout& block : layout.blocks)
    {
        const Matz& lhs = output.at(block.spin).at(block.kpoint);
        const Matz& rhs = input.at(block.spin).at(block.kpoint);
        for (int row = 0; row < block.rows; ++row)
        {
            for (int column = 0; column < block.columns; ++column)
            {
                const double magnitude = std::abs(
                    lhs(row, column) - rhs(row, column));
                square_sum += magnitude * magnitude;
                maximum = std::max(maximum, magnitude);
            }
        }
    }
    return {std::sqrt(square_sum), maximum};
}

SpinKHamiltonianMixer::SpinKHamiltonianMixer(MixingOptions options)
    : mixer_(std::move(options))
{
}

void SpinKHamiltonianMixer::initialize(const SpinKMatrixMap& grid_input)
{
    const MapLayout grid_layout = build_layout(grid_input, "grid");
    HamiltonianMixer next_mixer = mixer_;
    next_mixer.initialize(pack_map(grid_input, grid_layout));

    mixer_ = std::move(next_mixer);
    grid_layout_ = grid_layout;
    band_layout_.reset();
    initialized_ = true;
}

void SpinKHamiltonianMixer::initialize(const SpinKMatrixMap& grid_input,
                                       const SpinKMatrixMap& band_input)
{
    const MapLayout grid_layout = build_layout(grid_input, "grid");
    const MapLayout band_layout = build_layout(band_input, "band");
    HamiltonianMixer next_mixer = mixer_;
    next_mixer.initialize(pack_map(grid_input, grid_layout),
                          pack_map(band_input, band_layout));

    mixer_ = std::move(next_mixer);
    grid_layout_ = grid_layout;
    band_layout_ = band_layout;
    initialized_ = true;
}

SpinKHamiltonianMixResult SpinKHamiltonianMixer::mix(
    const SpinKMatrixMap& grid_output)
{
    return mix_impl(grid_output, std::nullopt);
}

SpinKHamiltonianMixResult SpinKHamiltonianMixer::mix(
    const SpinKMatrixMap& grid_output,
    const SpinKMatrixMap& band_output)
{
    return mix_impl(grid_output, band_output);
}

SpinKHamiltonianMixResult SpinKHamiltonianMixer::mix_impl(
    const SpinKMatrixMap& grid_output,
    const std::optional<SpinKMatrixMap>& band_output)
{
    if (!initialized_)
    {
        throw std::logic_error(
            "QSGW spin/k Hamiltonian mixer must be initialized before use");
    }
    require_matching_layout(grid_output, grid_layout_, "grid");
    if (band_layout_.has_value() != band_output.has_value())
    {
        throw std::invalid_argument(
            "QSGW band Hamiltonian presence changed during mixing");
    }
    if (band_output.has_value())
    {
        require_matching_layout(*band_output, *band_layout_, "band");
    }

    HamiltonianMixer next_mixer = mixer_;
    HamiltonianMixResult packed_result;
    if (band_output.has_value())
    {
        packed_result = next_mixer.mix(
            pack_map(grid_output, grid_layout_),
            pack_map(*band_output, *band_layout_));
    }
    else
    {
        packed_result = next_mixer.mix(pack_map(grid_output, grid_layout_));
    }

    SpinKHamiltonianMixResult result;
    result.grid = unpack_map(packed_result.grid, grid_layout_);
    if (packed_result.band.has_value())
    {
        result.band = unpack_map(*packed_result.band, *band_layout_);
    }
    result.decision = std::move(packed_result.decision);
    result.residual_l2 = packed_result.residual_l2;
    result.residual_max = packed_result.residual_max;

    mixer_ = std::move(next_mixer);
    return result;
}

} // namespace qsgw
} // namespace librpa_int
