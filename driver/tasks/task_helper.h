#pragma once

#include "../../src/core/gw_output.h"
#include "librpa.hpp"

using librpa_int::occupation_number_from_kpoint_weight;
using librpa_int::print_band_analysis;
using librpa_int::write_energy_qp;
using librpa_int::write_gw_band_spin_table;

//! Independent output-only copy retrieved through the public API; no WFC data.
librpa_int::MeanField meanfield_for_output(const librpa::MeanfieldData &snapshot);
