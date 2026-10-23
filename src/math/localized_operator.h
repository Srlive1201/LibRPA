#pragma once

#include <array>
#include <map>
#include <vector>

#include "matrix3.h"
#include "matrix_m.h"
#include "vector3_order.h"

namespace librpa_int
{

// Dense blocks with independent lattice-translation supports. Storage order,
// physical units and Fourier-completion conventions belong to the caller.
using LocalizedMatrix = std::map<Vector3_Order<int>, Matz>;

struct LocalizedAoOperators
{
    int dimension = 0;
    LocalizedMatrix overlap;
    std::array<LocalizedMatrix, 3> position;
    // Optional immutable reference H(R), one model per spin channel. Preserve
    // its support so callers can reconstruct only a current-state correction.
    std::vector<LocalizedMatrix> native_dft_hamiltonian;
};

// Shared localized-operator transforms. Fourier phases use fractional k and
// exp(+i 2*pi*k.R); Cartesian derivatives use row lattice vectors in Bohr.
namespace localized_operator
{
Matz zero_matrix(int n);
double maximum_entry(const Matz& matrix);
cplxdb phase(const Vector3_Order<double>& k, const Vector3_Order<int>& r);
Matz fourier(const LocalizedMatrix& blocks, const Vector3_Order<double>& k, int n,
             const Matrix3& lattice, int derivative = -1);
// Own the matrix before completing it when the input must remain unchanged:
// Matz ordinary copies share their storage.
Matz complete_hermitian_upper(Matz upper);
}  // namespace localized_operator

}  // namespace librpa_int
