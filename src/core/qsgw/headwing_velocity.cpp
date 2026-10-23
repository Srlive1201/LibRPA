#include "headwing_velocity.h"

#include <algorithm>
#include <cmath>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>

#include "../../utils/constants.h"

namespace librpa_int::qsgw
{
using localized_operator::complete_hermitian_upper;
using localized_operator::fourier;
using localized_operator::maximum_entry;
using localized_operator::phase;
using localized_operator::zero_matrix;

namespace
{
void diagonalize_localized_model(const Matz& hamiltonian, const Matz& overlap_matrix,
                                 std::vector<double>& eigenvalues, Matz& coefficients)
{
    const int n = hamiltonian.nr();
    coefficients = hamiltonian.copy();
    coefficients.swap_to_col_major();
    Matz overlap = overlap_matrix.copy();
    overlap.swap_to_col_major();
    eigenvalues.resize(n);
    const int itype = 1;
    const char jobz = 'V', uplo = 'U';
    int info = 0;
    int lwork = -1;
    cplxdb query;
    std::vector<double> real_work(std::max(1, 3 * n - 2));
    // Use the existing LAPACK ABI on column-major copies. The older
    // ComplexMatrix zhegv connector passes info by value, losing failure status.
    zhegv_(&itype, &jobz, &uplo, &n, coefficients.ptr(), &n, overlap.ptr(), &n, eigenvalues.data(),
           &query, &lwork, real_work.data(), &info);
    if (info != 0)
        throw std::runtime_error("Localized generalized eigensolve workspace query failed");
    lwork = std::max(std::max(1, 2 * n - 1), static_cast<int>(std::ceil(query.real())));
    std::vector<cplxdb> work(lwork);
    zhegv_(&itype, &jobz, &uplo, &n, coefficients.ptr(), &n, overlap.ptr(), &n, eigenvalues.data(),
           work.data(), &lwork, real_work.data(), &info);
    if (info != 0)
    {
        std::ostringstream message;
        message << "Localized generalized eigensolve failed (LAPACK info=" << info << ")";
        if (info > n) message << ": overlap is not positive definite";
        throw std::runtime_error(message.str());
    }
    coefficients.swap_to_row_major();
}

Matz ket_coefficients(const MeanField& mf, const int spin, const int ik)
{
    const int n = mf.get_n_bands();
    Matz u(n, n, MAJOR::ROW);
    if (mf.get_n_aos() * mf.get_n_spinor() != n)
        throw std::invalid_argument("QSGW live velocity needs the complete square AO reference");
    for (int spinor = 0; spinor < mf.get_n_spinor(); ++spinor)
    {
        const auto* block = mf.find_wfc(spin, spinor, ik);
        if (!block || block->nr != n || block->nc != mf.get_n_aos())
            throw std::invalid_argument(
                "QSGW live velocity requires replicated complete wavefunctions");
        for (int band = 0; band < n; ++band)
            for (int ao = 0; ao < mf.get_n_aos(); ++ao)
                u(ao * mf.get_n_spinor() + spinor, band) = (*block)(band, ao);
    }
    maximum_entry(u);
    return u;
}

double metric_error(const Matz& u, const Matz& s)
{
    Matz error = transpose(u, true) * s * u;
    for (int i = 0; i < error.nr(); ++i) error(i, i) -= 1.0;
    return maximum_entry(error);
}

double localized_hermiticity_error(const LocalizedMatrix& blocks, const int n)
{
    double result = 0.0;
    for (const auto& [r, block] : blocks)
    {
        if (block.nr() != n || block.nc() != n)
            throw std::invalid_argument("QSGW native DFT Hamiltonian has inconsistent dimensions");
        maximum_entry(block);
        const auto opposite = blocks.find(Vector3_Order<int>(-r.x, -r.y, -r.z));
        if (opposite == blocks.end())
            result = std::max(result, maximum_entry(block));
        else
        {
            if (opposite->second.nr() != n || opposite->second.nc() != n)
                throw std::invalid_argument(
                    "QSGW native DFT Hamiltonian has inconsistent dimensions");
            maximum_entry(opposite->second);
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < n; ++j)
                    result =
                        std::max(result, std::abs(block(i, j) - std::conj(opposite->second(j, i))));
        }
    }
    return result;
}

std::vector<int> images(const int index, const int size)
{
    if (size % 2 == 0 && index == size / 2) return {-size / 2, size / 2};
    return {index <= size / 2 ? index : index - size};
}
}  // namespace

LocalizedVelocityResult solve_localized_velocity(const LocalizedMatrix& hamiltonian,
                                                 const HeadwingAoOperators& operators,
                                                 const Vector3_Order<double>& kpoint_fractional,
                                                 const Matrix3& lattice_bohr)
{
    const int n = operators.dimension;
    if (n <= 0 || hamiltonian.empty() || operators.overlap.empty())
        throw std::invalid_argument("Localized velocity requires a nonempty square H/S model");
    const auto& k = kpoint_fractional;
    const auto& lattice = lattice_bohr;
    const std::array<double, 12> coordinates{k.x,         k.y,         k.z,         lattice.e11,
                                             lattice.e12, lattice.e13, lattice.e21, lattice.e22,
                                             lattice.e23, lattice.e31, lattice.e32, lattice.e33};
    for (const double value : coordinates)
        if (!std::isfinite(value))
            throw std::invalid_argument("Localized velocity requires finite k and lattice vectors");
    const Matz h = fourier(hamiltonian, k, n, lattice);
    const Matz s = complete_hermitian_upper(fourier(operators.overlap, k, n, lattice));
    std::vector<double> eigenvalues;
    Matz coefficients;
    diagonalize_localized_model(h, s, eigenvalues, coefficients);
    return compute_localized_velocity_in_frame(hamiltonian, operators, k, lattice, eigenvalues,
                                               coefficients);
}

HamiltonianReconstructionResult reconstruct_qsgw_hamiltonian(
    const MeanField& reference, const SpinKMatrixMap& current_hamiltonian,
    const HeadwingAoOperators& operators, const std::vector<Vector3_Order<double>>& kpoints,
    const Vector3_Order<int>& period, const Matrix3& lattice)
{
    const int n = reference.get_n_bands();
    const int nk = reference.get_n_kpoints();
    if (period.x <= 0 || period.y <= 0 || period.z <= 0 || nk != period.x * period.y * period.z ||
        static_cast<int>(kpoints.size()) != nk)
        throw std::invalid_argument(
            "QSGW Hamiltonian reconstruction requires a complete unshifted BvK mesh");
    if (operators.dimension != n || operators.overlap.empty())
        throw std::invalid_argument(
            "QSGW Hamiltonian reconstruction has inconsistent AO/reference dimensions");
    const bool preserve_native = !operators.native_dft_hamiltonian.empty();
    if (preserve_native &&
        static_cast<int>(operators.native_dft_hamiltonian.size()) != reference.get_n_spins())
        throw std::invalid_argument("QSGW native DFT Hamiltonian spin-channel count differs");
    std::set<std::array<int, 3>> points;
    for (const auto& k : kpoints)
    {
        std::array<int, 3> index{};
        const std::array<double, 3> fractional{k.x, k.y, k.z};
        const std::array<int, 3> sizes{period.x, period.y, period.z};
        for (int a = 0; a < 3; ++a)
        {
            const double x = fractional[a] * sizes[a];
            if (std::abs(x - std::round(x)) > 1e-6)
                throw std::invalid_argument(
                    "QSGW Hamiltonian reconstruction mesh is shifted or nonuniform");
            index[a] = (static_cast<int>(std::llround(x)) % sizes[a] + sizes[a]) % sizes[a];
        }
        if (!points.insert(index).second)
            throw std::invalid_argument(
                "QSGW Hamiltonian reconstruction mesh contains duplicate points");
    }

    std::vector<Matz> overlap;
    for (const auto& k : kpoints)
        overlap.push_back(complete_hermitian_upper(fourier(operators.overlap, k, n, lattice)));
    HamiltonianReconstructionResult result;
    result.hamiltonian.resize(reference.get_n_spins());
    for (int spin = 0; spin < reference.get_n_spins(); ++spin)
    {
        std::vector<Matz> h;
        std::vector<Matz> correction;
        auto& hr = result.hamiltonian[spin];
        if (preserve_native)
        {
            const auto& native = operators.native_dft_hamiltonian[spin];
            if (native.empty())
                throw std::invalid_argument("QSGW native DFT Hamiltonian has an empty spin model");
            result.native_dft_hermiticity_residual_ha = std::max(
                result.native_dft_hermiticity_residual_ha, localized_hermiticity_error(native, n));
            if (result.native_dft_hermiticity_residual_ha > 1e-9)
                throw std::runtime_error("QSGW native DFT Hamiltonian violates H(R)=H(-R)^dagger");
            // Matz has shared ordinary-copy storage: own every output block.
            for (const auto& [r, block] : native)
            {
                Matz owned = block.copy();
                owned.swap_to_row_major();
                hr.emplace(r, std::move(owned));
            }
        }
        for (int ik = 0; ik < nk; ++ik)
        {
            const Matz u = ket_coefficients(reference, spin, ik);
            result.reference_metric_residual =
                std::max(result.reference_metric_residual, metric_error(u, overlap[ik]));
            if (result.reference_metric_residual > 1e-6)
                throw std::runtime_error(
                    "QSGW Hamiltonian reconstruction reference metric exceeds tolerance");
            const Matz su = overlap[ik] * u;
            // The exported C need only satisfy the metric gate, not C^dagger S C = I
            // exactly. S C h C^dagger S would project back as M h M, amplifying
            // small export errors by high eigenvalues. Use the dual frame so
            // C^dagger H_AO C = h even at finite input precision. Invert the
            // near-identity metric instead of the potentially ill-conditioned C.
            const Matz metric = transpose(u, true) * su;
            Matz metric_inverse = metric.copy();
            metric_inverse.swap_to_row_major();
            std::vector<int> pivots(n);
            std::vector<cplxdb> work(n);
            int info = 0;
            LapackConnector::getrf(n, n, metric_inverse.ptr(), n, pivots.data(), info);
            if (info != 0) throw std::runtime_error("QSGW reference metric factorization failed");
            LapackConnector::getri(n, metric_inverse.ptr(), n, pivots.data(), work.data(), n, info);
            if (info != 0) throw std::runtime_error("QSGW reference metric inversion failed");
            const Matz dual = su * metric_inverse;
            h.push_back(dual * current_hamiltonian.at(spin).at(ik) * transpose(dual, true));
            if (preserve_native)
            {
                const Matz native =
                    fourier(operators.native_dft_hamiltonian[spin], kpoints[ik], n, lattice);
                Matz frame_error = transpose(u, true) * native * u;
                for (int band = 0; band < n; ++band)
                    frame_error(band, band) -= reference.get_eigenvals().at(spin)(ik, band);
                result.native_dft_reference_residual_ha =
                    std::max(result.native_dft_reference_residual_ha, maximum_entry(frame_error));
                correction.push_back(h.back() - native);
            }
        }
        if (result.native_dft_reference_residual_ha > 1e-4 / HA2EV)
            throw std::runtime_error(
                "QSGW native DFT Hamiltonian does not match the immutable reference on the mesh");
        const auto& samples = preserve_native ? correction : h;
        for (int ix = 0; ix < period.x; ++ix)
            for (int iy = 0; iy < period.y; ++iy)
                for (int iz = 0; iz < period.z; ++iz)
                {
                    const auto xs = images(ix, period.x), ys = images(iy, period.y),
                               zs = images(iz, period.z);
                    Matz block = zero_matrix(n);
                    const Vector3_Order<int> residue(ix, iy, iz);
                    for (int ik = 0; ik < nk; ++ik)
                    {
                        const cplxdb weight =
                            std::conj(phase(kpoints[ik], residue)) / static_cast<double>(nk);
                        for (std::size_t j = 0; j < block.size(); ++j)
                            block.ptr()[j] += weight * samples[ik].ptr()[j];
                    }
                    for (const int x : xs)
                        for (const int y : ys)
                            for (const int z : zs)
                            {
                                Matz image = block.copy();
                                const double weight = 1.0 / (xs.size() * ys.size() * zs.size());
                                for (std::size_t j = 0; j < image.size(); ++j)
                                    image.ptr()[j] *= weight;
                                const Vector3_Order<int> r(x, y, z);
                                const auto found = hr.find(r);
                                if (found == hr.end())
                                    hr.emplace(r, std::move(image));
                                else
                                    found->second += image;
                            }
                }
        for (int ik = 0; ik < nk; ++ik)
        {
            const Matz reconstructed = fourier(hr, kpoints[ik], n, lattice);
            result.fourier_roundtrip_residual_ha = std::max(result.fourier_roundtrip_residual_ha,
                                                            maximum_entry(reconstructed - h[ik]));
        }
    }
    if (result.reference_metric_residual > 1e-6 || result.fourier_roundtrip_residual_ha > 1e-9)
    {
        std::ostringstream message;
        message << std::scientific
                << "QSGW Hamiltonian reconstruction check failed: reference_metric="
                << result.reference_metric_residual
                << " FT_Ha=" << result.fourier_roundtrip_residual_ha;
        throw std::runtime_error(message.str());
    }
    return result;
}

LiveVelocityResult compute_live_headwing_velocity(const MeanField& reference, const MeanField& live,
                                                  const SpinKMatrixMap& current_hamiltonian,
                                                  const HeadwingAoOperators& operators,
                                                  const std::vector<Vector3_Order<double>>& kpoints,
                                                  const Vector3_Order<int>& period,
                                                  const Matrix3& lattice)
{
    const int n = reference.get_n_bands();
    const int nk = reference.get_n_kpoints();
    if (live.get_n_bands() != n || live.get_n_kpoints() != nk ||
        live.get_n_spins() != reference.get_n_spins() ||
        live.get_n_aos() != reference.get_n_aos() ||
        live.get_n_spinor() != reference.get_n_spinor())
        throw std::invalid_argument(
            "QSGW current-HR velocity has inconsistent AO/mean-field dimensions");
    auto reconstruction = reconstruct_qsgw_hamiltonian(reference, current_hamiltonian, operators,
                                                       kpoints, period, lattice);
    LiveVelocityResult result;
    result.hamiltonian = std::move(reconstruction.hamiltonian);
    result.reference_metric_residual = reconstruction.reference_metric_residual;
    result.fourier_roundtrip_residual_ha = reconstruction.fourier_roundtrip_residual_ha;
    result.native_dft_reference_residual_ha = reconstruction.native_dft_reference_residual_ha;
    result.native_dft_hermiticity_residual_ha = reconstruction.native_dft_hermiticity_residual_ha;
    result.velocity.resize(reference.get_n_spins());
    for (int spin = 0; spin < reference.get_n_spins(); ++spin)
    {
        const auto& hr = result.hamiltonian[spin];
        result.velocity[spin].resize(nk);
        for (int ik = 0; ik < nk; ++ik)
        {
            const Matz u = ket_coefficients(live, spin, ik);
            std::vector<double> eigenvalues(n);
            for (int i = 0; i < n; ++i) eigenvalues[i] = live.get_eigenvals().at(spin)(ik, i);
            const auto projected = compute_localized_velocity_in_frame(hr, operators, kpoints[ik],
                                                                       lattice, eigenvalues, u);
            result.live_metric_residual =
                std::max(result.live_metric_residual, projected.metric_residual);
            result.hamiltonian_frame_residual_ha = std::max(
                result.hamiltonian_frame_residual_ha, projected.hamiltonian_frame_residual_ha);
            result.velocity_upper_diagonal_imaginary_ha_bohr =
                std::max(result.velocity_upper_diagonal_imaginary_ha_bohr,
                         projected.velocity_upper_diagonal_imaginary_ha_bohr);
            for (int a = 0; a < 3; ++a)
            {
                ComplexMatrix v(n, n);
                for (int i = 0; i < n; ++i)
                    for (int j = 0; j < n; ++j) v(i, j) = projected.velocity_ha_bohr[a](i, j);
                result.velocity[spin][ik].push_back(std::move(v));
            }
        }
    }
    if (result.reference_metric_residual > 1e-6 || result.live_metric_residual > 1e-6 ||
        result.hamiltonian_frame_residual_ha > 1e-4 / HA2EV ||
        result.fourier_roundtrip_residual_ha > 1e-9 ||
        result.velocity_upper_diagonal_imaginary_ha_bohr > 1e-8)
    {
        std::ostringstream message;
        message << std::scientific << "QSGW current-HR velocity check failed: reference_metric="
                << result.reference_metric_residual
                << " live_metric=" << result.live_metric_residual
                << " frame_Ha=" << result.hamiltonian_frame_residual_ha
                << " FT_Ha=" << result.fourier_roundtrip_residual_ha
                << " diagonal_imag_HaBohr=" << result.velocity_upper_diagonal_imaginary_ha_bohr;
        throw std::runtime_error(message.str());
    }
    return result;
}
}  // namespace librpa_int::qsgw
