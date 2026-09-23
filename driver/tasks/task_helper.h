#pragma once

#include <vector>

#include "../../src/core/meanfield.h"
#include "../../src/utils/base_utility.h"

void print_band_analysis(const librpa_int::MeanField &mf,
                         const std::vector<librpa_int::Vector3_Order<double>> &kfrac_output,
                         const std::vector<int> &output_to_input_kpoint,
                         const std::vector<librpa_int::matrix> &vxc,
                         const std::vector<double> &vexx,
                         const std::vector<librpa_int::cplxdb> &sigc,
                         int i_state_low,
                         int n_states_calc);

void write_energy_qp(const librpa_int::MeanField &mf,
                     const std::vector<librpa_int::Vector3_Order<double>> &kfrac_output,
                     const std::vector<int> &output_to_input_kpoint,
                     const std::vector<librpa_int::matrix> &vxc,
                     const std::vector<double> &vexx,
                     const std::vector<librpa_int::cplxdb> &sigc,
                     int n_kpoints_data,
                     int i_state_low,
                     int n_states_calc,
                     double occupation_scale);
