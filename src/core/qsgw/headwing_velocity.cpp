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
namespace
{
constexpr double two_pi = 6.283185307179586476925286766559;

Matz zero_matrix(const int n)
{
    Matz result(n, n, MAJOR::ROW);
    std::fill(result.ptr(), result.ptr() + result.size(), cplxdb{});
    return result;
}

double maximum_entry(const Matz& matrix)
{
    double result = 0.0;
    for (std::size_t i = 0; i < matrix.size(); ++i)
    {
        if (!std::isfinite(matrix.ptr()[i].real()) || !std::isfinite(matrix.ptr()[i].imag()))
            throw std::runtime_error("QSGW live velocity contains non-finite matrix data");
        result = std::max(result, std::abs(matrix.ptr()[i]));
    }
    return result;
}

cplxdb phase(const Vector3_Order<double>& k, const Vector3_Order<int>& r)
{
    return std::polar(1.0, two_pi * (k.x * r.x + k.y * r.y + k.z * r.z));
}

Matz fourier(const LocalizedMatrix& blocks, const Vector3_Order<double>& k, const int n,
             const Matrix3& lattice, const int derivative = -1)
{
    Matz result = zero_matrix(n);
    for (const auto& [r, matrix] : blocks)
    {
        if (matrix.nr() != n || matrix.nc() != n)
            throw std::invalid_argument("QSGW live velocity localized basis dimensions differ");
        cplxdb weight = phase(k, r);
        if (derivative >= 0)
        {
            const std::array<double, 3> distance{
                r.x * lattice.e11 + r.y * lattice.e21 + r.z * lattice.e31,
                r.x * lattice.e12 + r.y * lattice.e22 + r.z * lattice.e32,
                r.x * lattice.e13 + r.y * lattice.e23 + r.z * lattice.e33};
            weight *= cplxdb(0.0, distance.at(derivative));
        }
        maximum_entry(matrix);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) result(i, j) += weight * matrix(i, j);
    }
    return result;
}

// S/r use the audited PyATB upper-triangle storage and completion convention.
Matz complete_hermitian_upper(Matz upper)
{
    for (int i = 0; i < upper.nr(); ++i)
    {
        upper(i, i) = upper(i, i).real();
        for (int j = i + 1; j < upper.nc(); ++j) upper(j, i) = std::conj(upper(i, j));
    }
    return upper;
}

struct FourierVelocityModel
{
    Matz hamiltonian;
    Matz overlap;
    std::array<Matz, 3> hamiltonian_derivative;
    std::array<Matz, 3> overlap_derivative;
    std::array<Matz, 3> position;
};

FourierVelocityModel fourier_velocity_model(const LocalizedMatrix& hamiltonian,
                                            const HeadwingAoOperators& operators,
                                            const Vector3_Order<double>& k, const Matrix3& lattice)
{
    const int n = operators.dimension;
    if (n <= 0 || hamiltonian.empty() || operators.overlap.empty())
        throw std::invalid_argument("Localized velocity requires a nonempty square H/S model");
    const std::array<double, 12> coordinates{k.x,         k.y,         k.z,         lattice.e11,
                                             lattice.e12, lattice.e13, lattice.e21, lattice.e22,
                                             lattice.e23, lattice.e31, lattice.e32, lattice.e33};
    for (const double value : coordinates)
        if (!std::isfinite(value))
            throw std::invalid_argument("Localized velocity requires finite k and lattice vectors");
    FourierVelocityModel result;
    result.hamiltonian = fourier(hamiltonian, k, n, lattice);
    result.overlap = complete_hermitian_upper(fourier(operators.overlap, k, n, lattice));
    for (int a = 0; a < 3; ++a)
    {
        result.hamiltonian_derivative[a] = fourier(hamiltonian, k, n, lattice, a);
        result.overlap_derivative[a] =
            complete_hermitian_upper(fourier(operators.overlap, k, n, lattice, a));
        Matz r = fourier(operators.position[a], k, n, lattice);
        const Matz raw_r = r.copy();
        for (int i = 0; i < n; ++i)
            for (int j = 0; j <= i; ++j)
                r(i, j) =
                    std::conj(raw_r(j, i) + cplxdb(0.0, 1.0) * result.overlap_derivative[a](j, i));
        result.position[a] = std::move(r);
    }
    return result;
}

LocalizedVelocityResult project_velocity(const FourierVelocityModel& model,
                                         const std::vector<double>& eigenvalues,
                                         const Matz& coefficients)
{
    const int n = model.hamiltonian.nr();
    if (static_cast<int>(eigenvalues.size()) != n || coefficients.nr() != n ||
        coefficients.nc() != n)
        throw std::invalid_argument("Localized velocity requires complete matching E/C dimensions");
    for (const double energy : eigenvalues)
        if (!std::isfinite(energy))
            throw std::invalid_argument("Localized velocity contains a non-finite eigenvalue");
    maximum_entry(coefficients);
    LocalizedVelocityResult result;
    result.eigenvalues_ha = eigenvalues;
    result.ket_coefficients = coefficients.copy();
    result.ket_coefficients.swap_to_row_major();
    result.overlap = model.overlap.copy();
    const Matz& u = result.ket_coefficients;
    const Matz ud = transpose(u, true);
    Matz metric = ud * model.overlap * u;
    for (int i = 0; i < n; ++i) metric(i, i) -= 1.0;
    result.metric_residual = maximum_entry(metric);
    Matz state_h = ud * model.hamiltonian * u;
    for (int i = 0; i < n; ++i) state_h(i, i) -= eigenvalues[i];
    result.hamiltonian_frame_residual_ha = maximum_entry(state_h);
    result.hamiltonian_hermiticity_residual_ha =
        maximum_entry(model.hamiltonian - transpose(model.hamiltonian, true));
    for (int a = 0; a < 3; ++a)
    {
        const Matz dh = ud * model.hamiltonian_derivative[a] * u;
        const Matz ds = ud * model.overlap_derivative[a] * u;
        const Matz r = ud * model.position[a] * u;
        Matz v = zero_matrix(n);
        for (int i = 0; i < n; ++i)
            for (int j = i; j < n; ++j)
            {
                const cplxdb value = dh(i, j) - eigenvalues[i] * ds(i, j) +
                                     cplxdb(0.0, eigenvalues[i] - eigenvalues[j]) * r(i, j);
                if (!std::isfinite(value.real()) || !std::isfinite(value.imag()))
                    throw std::runtime_error("Localized velocity contains non-finite values");
                if (i == j)
                    result.velocity_upper_diagonal_imaginary_ha_bohr = std::max(
                        result.velocity_upper_diagonal_imaginary_ha_bohr, std::abs(value.imag()));
                v(i, j) = i == j ? cplxdb(value.real(), 0.0) : value;
                v(j, i) = std::conj(v(i, j));
            }
        result.velocity_ha_bohr[a] = std::move(v);
    }
    // Retain the existing live-path gates. The Hermiticity residual is reported
    // independently; an inconsistent full H also enters the eigenframe check.
    if (result.metric_residual > 1e-6 || result.hamiltonian_frame_residual_ha > 1e-4 / HA2EV ||
        result.velocity_upper_diagonal_imaginary_ha_bohr > 1e-8)
    {
        std::ostringstream message;
        message << std::scientific
                << "Localized velocity check failed: metric=" << result.metric_residual
                << " frame_Ha=" << result.hamiltonian_frame_residual_ha
                << " diagonal_imag_HaBohr=" << result.velocity_upper_diagonal_imaginary_ha_bohr;
        throw std::runtime_error(message.str());
    }
    return result;
}

void diagonalize_localized_model(const FourierVelocityModel& model,
                                 std::vector<double>& eigenvalues, Matz& coefficients)
{
    const int n = model.hamiltonian.nr();
    coefficients = model.hamiltonian.copy();
    coefficients.swap_to_col_major();
    Matz overlap = model.overlap.copy();
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
    const auto model =
        fourier_velocity_model(hamiltonian, operators, kpoint_fractional, lattice_bohr);
    std::vector<double> eigenvalues;
    Matz coefficients;
    diagonalize_localized_model(model, eigenvalues, coefficients);
    return project_velocity(model, eigenvalues, coefficients);
}

LocalizedVelocityResult compute_localized_velocity_in_frame(
    const LocalizedMatrix& hamiltonian, const HeadwingAoOperators& operators,
    const Vector3_Order<double>& kpoint_fractional, const Matrix3& lattice_bohr,
    const std::vector<double>& eigenvalues_ha, const Matz& ket_coefficients)
{
    return project_velocity(
        fourier_velocity_model(hamiltonian, operators, kpoint_fractional, lattice_bohr),
        eigenvalues_ha, ket_coefficients);
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
