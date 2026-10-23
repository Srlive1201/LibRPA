#pragma once

#include <iosfwd>
#include <vector>

#include "../../math/vector3_order.h"
#include "../meanfield.h"
#include "matrix_map.h"

namespace librpa_int
{
namespace qsgw
{

// Write one spin channel using the historical QSGW/KS/EXX band-table layout.
// Matrix inputs and MeanField eigenvalues are in Hartree; table energies are eV.
void write_qsgw_band_spin_tables(std::ostream& ks_output, std::ostream& exx_output,
                                 std::ostream& qsgw_output, const MeanField& live_band,
                                 const MeanField& reference_band,
                                 const std::vector<Vector3_Order<double>>& kpoints,
                                 const SpinKMatrixMap& reference_hamiltonian,
                                 const SpinKMatrixMap& dft_vxc, const SpinKMatrixMap& exchange,
                                 int spin, double chemical_potential_ha);

}  // namespace qsgw
}  // namespace librpa_int

namespace librpa_int::qsgw
{
double occupation_number_from_kpoint_weight(double weighted_occupation, double kpoint_weight);

void write_energy_qp(std::ostream &ofs, const librpa_int::MeanField &mf,
                     const std::vector<librpa_int::Vector3_Order<double>> &kfrac_output,
                     const std::vector<int> &output_to_input_kpoint,
                     const std::vector<librpa_int::matrix> &vxc, const std::vector<double> &vexx,
                     const std::vector<librpa_int::cplxdb> &sigc, const int n_kpoints_data,
                     const int i_state_low, const int n_states_calc,
                     const std::vector<double> &kpoint_weights);

void write_gw_band_spin_table(std::ostream &output, const librpa_int::MeanField &mf,
                              const std::vector<librpa_int::Vector3_Order<double>> &kfrac,
                              const std::vector<librpa_int::matrix> &vxc,
                              const std::vector<double> &vexx,
                              const std::vector<librpa_int::cplxdb> &sigc, int spin,
                              int i_state_low, int n_states_calc);
}  // namespace librpa_int::qsgw
