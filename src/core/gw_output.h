#pragma once

#include <iosfwd>
#include <string>
#include <vector>

#include "meanfield.h"
#include "../utils/base_utility.h"

namespace librpa_int
{

void print_band_analysis(const librpa_int::MeanField &mf,
                         const std::vector<librpa_int::Vector3_Order<double>> &kfrac_output,
                         const std::vector<int> &output_to_input_kpoint,
                         const std::vector<librpa_int::matrix> &vxc,
                         const std::vector<double> &vexx,
                         const std::vector<librpa_int::cplxdb> &sigc,
                         int i_state_low,
                         int n_states_calc);

double occupation_number_from_kpoint_weight(double weighted_occupation, double kpoint_weight);

void write_energy_qp(const librpa_int::MeanField &mf,
                     const std::vector<librpa_int::Vector3_Order<double>> &kfrac_output,
                     const std::vector<int> &output_to_input_kpoint,
                     const std::vector<librpa_int::matrix> &vxc, const std::vector<double> &vexx,
                     const std::vector<librpa_int::cplxdb> &sigc, int n_kpoints_data,
                     int i_state_low, int n_states_calc, const std::vector<double> &kpoint_weights,
                     const std::string &filename = "energy_qp");

// Native G0W0 path table: k index/coordinates, followed by occupation/energy
// pairs. Inputs are in Hartree; output energies are in eV, without a header.
void write_gw_band_spin_table(std::ostream &output, const librpa_int::MeanField &mf,
                              const std::vector<librpa_int::Vector3_Order<double>> &kfrac,
                              const std::vector<librpa_int::matrix> &vxc,
                              const std::vector<double> &vexx,
                              const std::vector<librpa_int::cplxdb> &sigc, int spin,
                              int i_state_low, int n_states_calc);

void write_energy_qp(std::ostream &ofs, const librpa_int::MeanField &mf,
                     const std::vector<librpa_int::Vector3_Order<double>> &kfrac_output,
                     const std::vector<int> &output_to_input_kpoint,
                     const std::vector<librpa_int::matrix> &vxc, const std::vector<double> &vexx,
                     const std::vector<librpa_int::cplxdb> &sigc, const int n_kpoints_data,
                     const int i_state_low, const int n_states_calc,
                     const std::vector<double> &kpoint_weights);

}  // namespace librpa_int
