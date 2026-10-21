#pragma once
#include "../../src/core/gw_output.h"
#include "librpa.hpp"
using librpa_int::occupation_number_from_kpoint_weight;
using librpa_int::print_band_analysis;
librpa_int::MeanField meanfield_for_output(const librpa::MeanfieldData&);
using librpa_int::write_energy_qp;
