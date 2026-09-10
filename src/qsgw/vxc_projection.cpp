#include "vxc_projection.h"
#include <cmath>
#include <stdexcept>

namespace librpa_int
{
namespace qsgw
{
namespace
{
constexpr double hermitian_tolerance = 1.0e-10;
bool finite_complex(const cplxdb value)
{
    return std::isfinite(value.real()) && std::isfinite(value.imag());
}

void validate_hermitian_matrix(const Matz& matrix,
                               const int expected_dimension,
                               const std::string& label)
{
    if (matrix.nr() != expected_dimension ||
        matrix.nc() != expected_dimension)
    {
        throw std::invalid_argument(label + " has an invalid shape");
    }
    for (int row = 0; row < matrix.nr(); ++row)
    {
        for (int column = 0; column < matrix.nc(); ++column)
        {
            if (!finite_complex(matrix(row, column)))
            {
                throw std::invalid_argument(label + " contains non-finite data");
            }
            if (std::abs(matrix(row, column) -
                         std::conj(matrix(column, row))) >
                hermitian_tolerance)
            {
                throw std::invalid_argument(label + " is not Hermitian");
            }
        }
    }
}

void validate_reference_wfc(const MeanField& reference,
                            const int spin,
                            const int kpoint)
{
    if (!reference.initialized() || spin < 0 ||
        spin >= reference.get_n_spins() || kpoint < 0 ||
        kpoint >= reference.get_n_kpoints())
    {
        throw std::invalid_argument(
            "QSGW Vxc projection received an invalid mean-field index");
    }
    for (int spinor = 0; spinor < reference.get_n_spinor(); ++spinor)
    {
        const ComplexMatrix* wfc = reference.find_wfc(spin, spinor, kpoint);
        if (wfc == nullptr || wfc->nr != reference.get_n_bands() ||
            wfc->nc != reference.get_n_aos())
        {
            throw std::invalid_argument(
                "QSGW Vxc projection reference wavefunction is incomplete");
        }
        for (int index = 0; index < wfc->size; ++index)
        {
            if (!finite_complex(wfc->c[index]))
            {
                throw std::invalid_argument(
                    "QSGW Vxc projection wavefunction contains non-finite data");
            }
        }
    }
}

}
Matz project_vxc_nao_to_fixed_basis(const Matz& vxc_nao,
                                    const MeanField& reference,
                                    const int spin,
                                    const int kpoint)
{
    validate_reference_wfc(reference, spin, kpoint);
    validate_hermitian_matrix(vxc_nao, reference.get_n_aos(),
                              "QSGW NAO Vxc matrix");
    Matz result(reference.get_n_bands(), reference.get_n_bands(),
                MAJOR::ROW);
    for (int bra = 0; bra < reference.get_n_bands(); ++bra)
    {
        for (int ket = 0; ket < reference.get_n_bands(); ++ket)
        {
            for (int spinor = 0; spinor < reference.get_n_spinor(); ++spinor)
            {
                const ComplexMatrix& wfc =
                    reference.get_eigenvectors()
                        .at(spin)
                        .at(spinor)
                        .at(kpoint);
                for (int row = 0; row < reference.get_n_aos(); ++row)
                {
                    for (int column = 0; column < reference.get_n_aos();
                         ++column)
                    {
                        result(bra, ket) +=
                            std::conj(wfc(bra, row)) * vxc_nao(row, column) *
                            wfc(ket, column);
                    }
                }
            }
        }
    }
    validate_hermitian_matrix(result, reference.get_n_bands(),
                              "QSGW projected Vxc matrix");
    return result;
}

Matz prepare_vxc_in_fixed_state_basis(const Matz& input,
                                      const VxcBasis basis,
                                      const MeanField& reference,
                                      const int spin,
                                      const int kpoint)
{
    validate_reference_wfc(reference, spin, kpoint);
    if (basis == VxcBasis::Nao)
    {
        return project_vxc_nao_to_fixed_basis(
            input, reference, spin, kpoint);
    }
    if (basis == VxcBasis::State)
    {
        validate_hermitian_matrix(input, reference.get_n_bands(),
                                  "QSGW state-basis Vxc matrix");
        return input.copy();
    }
    throw std::invalid_argument("QSGW Vxc basis is invalid");
}

} // namespace qsgw
} // namespace librpa_int
