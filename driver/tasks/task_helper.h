#pragma once

#include <complex>
#include <vector>

#include "../../src/math/matrix.h"
#include "../../src/math/vector3_order.h"
#include "../../src/utils/base_utility.h"
#include "librpa.hpp"

//! Borrowed public output view with owning storage for layout conversions only.
struct GwOutputView
{
    std::vector<double> fractional_k, vxc;
    LibrpaGwOutputInput input{};

    GwOutputView(const librpa::MeanfieldData &snapshot, bool band_path,
                 const std::vector<librpa_int::Vector3_Order<double>> &kfrac,
                 const std::vector<int> &mapping, const std::vector<librpa_int::matrix> &vxc,
                 const std::vector<double> &vexx, const std::vector<std::complex<double>> &sigc,
                 int state_low, int state_high, const std::vector<double> &weights = {});
};

//! Flat index in a public owning spectrum/occupation snapshot.
std::size_t meanfield_output_index(const librpa::MeanfieldData &snapshot, int spin, int k,
                                   int state);
