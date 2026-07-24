#pragma once

#include "matrix_map.h"

#include "../core/meanfield.h"
#include "../math/vector3_order.h"

#include <iosfwd>
#include <vector>

namespace librpa_int
{
namespace qsgw
{

// Write one spin channel using the historical QSGW/KS/EXX band-table layout.
// Matrix inputs and MeanField eigenvalues are in Hartree; table energies are eV.
void write_qsgw_band_spin_tables(
    std::ostream& ks_output,
    std::ostream& exx_output,
    std::ostream& qsgw_output,
    const MeanField& live_band,
    const MeanField& reference_band,
    const std::vector<Vector3_Order<double>>& kpoints,
    const SpinKMatrixMap& reference_hamiltonian,
    const SpinKMatrixMap& dft_vxc,
    const SpinKMatrixMap& exchange,
    int spin,
    double chemical_potential_ha);

} // namespace qsgw
} // namespace librpa_int
