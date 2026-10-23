#pragma once

#include <string>

#include "../../src/math/localized_operator.h"

namespace librpa::reader
{
// Native ABACUS compact/legacy/HContainer CSR: S dimensionless, r in Bohr.
// QSGW supplies its current H in memory. G0W0 and native-baseline QSGW export
// read the original DFT H below.
librpa_int::LocalizedAoOperators read_headwing_ao_operators(const std::string& overlap_file,
                                                            const std::string& position_file,
                                                            double position_scale_bohr = 1.0);
// Original DFT H(R): retain all complex entries and independent R support.
// The driver supplies the explicit conversion from input energy to Hartree.
librpa_int::LocalizedMatrix read_headwing_hamiltonian_csr(const std::string& filename,
                                                          double energy_scale_ha);
// Diagnostic model actually differentiated for one entering QSGW update.
// Full complex CSR values are in eV; no entry threshold is applied.
void write_headwing_hamiltonian_csr(const std::string& filename,
                                    const librpa_int::LocalizedMatrix& blocks);
// Complete mixed/post-cut QSGW operator: Ha input, HContainer complex CSR17 in Ry.
// Retain the complete H/S/r support union, padding missing H blocks with exact
// zeros without changing the original S/r operators. Publish by atomic rename.
// Collinear channels have spin_directions=2; scalar/SOC use 1. Index is one-based.
void write_qsgw_hamiltonian_csr(const std::string& filename,
                                const librpa_int::LocalizedMatrix& blocks,
                                const librpa_int::LocalizedAoOperators& operators, int iteration,
                                int spin_directions = 1, int spin_index = 1);
}  // namespace librpa::reader
