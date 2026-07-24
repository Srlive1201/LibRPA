#pragma once

#include "matrix_map.h"
#include "mixing.h"

#include <optional>
#include <vector>

namespace librpa_int
{
namespace qsgw
{

struct SpinKHamiltonianMixResult
{
    SpinKMatrixMap grid;
    std::optional<SpinKMatrixMap> band;
    MixingDecision decision;
    double residual_l2 = 0.0;
    double residual_max = 0.0;
};

struct SpinKHamiltonianResidual
{
    double l2 = 0.0;
    double maximum = 0.0;
};

SpinKHamiltonianResidual measure_spin_k_hamiltonian_residual(
    const SpinKMatrixMap& output,
    const SpinKMatrixMap& input);

namespace detail
{

struct HamiltonianMatrixBlockLayout
{
    int spin = 0;
    int kpoint = 0;
    int rows = 0;
    int columns = 0;
    MAJOR major = MAJOR::ROW;
    int packed_offset = 0;
};

struct HamiltonianMapLayout
{
    std::vector<HamiltonianMatrixBlockLayout> blocks;
    int packed_size = 0;
};

} // namespace detail

class SpinKHamiltonianMixer
{
public:
    explicit SpinKHamiltonianMixer(MixingOptions options = {});

    void initialize(const SpinKMatrixMap& grid_input);
    void initialize(const SpinKMatrixMap& grid_input,
                    const SpinKMatrixMap& band_input);

    SpinKHamiltonianMixResult mix(const SpinKMatrixMap& grid_output);
    SpinKHamiltonianMixResult mix(const SpinKMatrixMap& grid_output,
                                 const SpinKMatrixMap& band_output);

private:
    HamiltonianMixer mixer_;
    bool initialized_ = false;
    detail::HamiltonianMapLayout grid_layout_;
    std::optional<detail::HamiltonianMapLayout> band_layout_;

    SpinKHamiltonianMixResult mix_impl(
        const SpinKMatrixMap& grid_output,
        const std::optional<SpinKMatrixMap>& band_output);
};

} // namespace qsgw
} // namespace librpa_int
