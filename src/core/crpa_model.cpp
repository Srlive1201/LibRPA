#include "crpa_model.h"

#include <cmath>

#include "../utils/constants.h"
#include "../utils/error.h"

namespace librpa_int
{
namespace
{
void require_isometry(const ComplexMatrix& rotation)
{
    if (rotation.nr <= 0 || rotation.nc <= 0 || rotation.nc > rotation.nr)
        throw LIBRPA_INVALID_ARGUMENT("orbital rotation needs a nonempty isometry");
    const auto gram = transpose(rotation, true) * rotation;
    for (int i = 0; i < gram.nr; ++i)
        for (int j = 0; j < gram.nc; ++j)
            if (!std::isfinite(std::abs(gram(i, j))) ||
                std::abs(gram(i, j) - (i == j ? 1.0 : 0.0)) > 1e-10)
                throw LIBRPA_INVALID_ARGUMENT("orbital rotation is not an isometry");
}

void require_tensor_dimensions(const ComplexMatrix& tensor, int n_orbitals)
{
    if (n_orbitals <= 0 || tensor.nr != n_orbitals * n_orbitals || tensor.nc != tensor.nr)
        throw LIBRPA_INVALID_ARGUMENT("ordered cRPA tensor has inconsistent orbital dimensions");
}
}  // namespace

ComplexMatrix project_crpa_one_body(const ComplexMatrix& projector,
                                    const std::vector<double>& diagonal)
{
    if (projector.nc != diagonal.size())
        throw LIBRPA_INVALID_ARGUMENT("projected operator has inconsistent band dimensions");
    ComplexMatrix weighted(projector);
    for (int band = 0; band < projector.nc; ++band)
    {
        if (!std::isfinite(diagonal[band]))
            throw LIBRPA_INVALID_ARGUMENT("projected operator has nonfinite diagonal entries");
        for (int orbital = 0; orbital < projector.nr; ++orbital)
            weighted(orbital, band) *= diagonal[band];
    }
    return weighted * transpose(projector, true);
}

ComplexMatrix project_crpa_orbitals(const ComplexMatrix& coefficients,
                                    const ComplexMatrix& spatial_overlap, const ComplexMatrix& phi,
                                    int n_spinor)
{
    if ((n_spinor != 1 && n_spinor != 2) || spatial_overlap.nr != spatial_overlap.nc ||
        coefficients.nr != spatial_overlap.nr * n_spinor || phi.nr != coefficients.nr)
        throw LIBRPA_INVALID_ARGUMENT("orbital projection has inconsistent spinor dimensions");
    const auto overlap = crpa_spinor_overlap(spatial_overlap, n_spinor);
    const auto gram = transpose(phi, true) * overlap * phi;
    for (int i = 0; i < phi.nc; ++i)
        for (int j = 0; j < phi.nc; ++j)
            if (!std::isfinite(std::abs(gram(i, j))) ||
                std::abs(gram(i, j) - (i == j ? 1.0 : 0.0)) > 1e-10)
                throw LIBRPA_INVALID_ARGUMENT("explicit orbital frame is not metric-orthonormal");
    return transpose(phi, true) * overlap * coefficients;
}

ComplexMatrix rotate_crpa_projector(const ComplexMatrix& projector, const ComplexMatrix& rotation)
{
    require_isometry(rotation);
    if (projector.nr != rotation.nr)
        throw LIBRPA_INVALID_ARGUMENT("projector rotation has inconsistent orbital dimensions");
    return transpose(rotation, true) * projector;
}

ComplexMatrix rotate_crpa_tensor(const ComplexMatrix& tensor, const ComplexMatrix& left_rotation,
                                 const ComplexMatrix& right_rotation)
{
    const auto pair_transform = [](const ComplexMatrix& rotation)
    {
        require_isometry(rotation);
        ComplexMatrix transform(rotation.nc * rotation.nc, rotation.nr * rotation.nr);
        for (int a = 0; a < rotation.nc; ++a)
            for (int b = 0; b < rotation.nc; ++b)
                for (int i = 0; i < rotation.nr; ++i)
                    for (int j = 0; j < rotation.nr; ++j)
                        transform(a * rotation.nc + b, i * rotation.nr + j) =
                            std::conj(rotation(i, a)) * rotation(j, b);
        return transform;
    };
    if (tensor.nr != left_rotation.nr * left_rotation.nr ||
        tensor.nc != right_rotation.nr * right_rotation.nr)
        throw LIBRPA_INVALID_ARGUMENT("interaction rotation has inconsistent dimensions");
    // The stored density tensor is bilinear in ordered densities, not a Hermitian pair matrix.
    return pair_transform(left_rotation) * tensor *
           transpose(pair_transform(right_rotation), false);
}

ComplexMatrix crpa_tensor_to_fermion(const ComplexMatrix& tensor, int n_orbitals)
{
    require_tensor_dimensions(tensor, n_orbitals);
    ComplexMatrix fermion(tensor.nr, tensor.nc);
    for (int i = 0; i < n_orbitals; ++i)
        for (int j = 0; j < n_orbitals; ++j)
            for (int k = 0; k < n_orbitals; ++k)
                for (int l = 0; l < n_orbitals; ++l)
                    fermion(i * n_orbitals + j, k * n_orbitals + l) =
                        tensor(i * n_orbitals + k, j * n_orbitals + l);
    return fermion;
}

std::array<std::complex<double>, 3> crpa_tensor_averages(const ComplexMatrix& tensor,
                                                         int n_orbitals)
{
    require_tensor_dimensions(tensor, n_orbitals);
    std::array<std::complex<double>, 3> averages{};
    for (int m = 0; m < n_orbitals; ++m)
    {
        averages[0] += tensor(m * n_orbitals + m, m * n_orbitals + m) / double(n_orbitals);
        for (int n = 0; n < n_orbitals; ++n)
            if (m != n)
            {
                const double count = n_orbitals * (n_orbitals - 1);
                averages[1] += tensor(m * n_orbitals + m, n * n_orbitals + n) / count;
                averages[2] += tensor(m * n_orbitals + n, n * n_orbitals + m) / count;
            }
    }
    return averages;
}

std::pair<ComplexMatrix, ComplexMatrix> average_crpa_local_model(
    const CrpaModelSnapshot& model, bool parent, int spin, int site,
    const std::vector<double>& k_weights)
{
    const auto& group = (parent ? model.parent_sites : model.output_sites).at(site);
    const auto& states = model.spin_k.at(spin);
    if (k_weights.size() != states.size())
        throw LIBRPA_INVALID_ARGUMENT("local cRPA model has inconsistent k weights");
    ComplexMatrix h(group.n_orbitals, group.n_orbitals), density(h);
    for (int k = 0; k < states.size(); ++k)
    {
        const auto& hk = parent ? states[k].parent_h : states[k].output_h;
        const auto& nk = parent ? states[k].parent_density : states[k].output_density;
        for (int a = 0; a < group.n_orbitals; ++a)
            for (int b = 0; b < group.n_orbitals; ++b)
            {
                h(a, b) += k_weights[k] * hk(group.orb_start + a, group.orb_start + b);
                density(a, b) += k_weights[k] * nk(group.orb_start + a, group.orb_start + b);
            }
    }
    return {std::move(h), std::move(density)};
}

ComplexMatrix fourier_crpa_model_hopping(const CrpaModelSnapshot& model, bool parent, int spin,
                                         const std::vector<Vector3_Order<double>>& k_fractional,
                                         const std::vector<double>& k_weights,
                                         const Vector3_Order<int>& lattice_R)
{
    const auto& states = model.spin_k.at(spin);
    if (states.empty() || k_fractional.size() != states.size() || k_weights.size() != states.size())
        throw LIBRPA_INVALID_ARGUMENT("cRPA hopping has inconsistent k sampling");
    const int n = parent ? states.front().parent_h.nr : states.front().output_h.nr;
    ComplexMatrix h(n, n);
    for (int k = 0; k < states.size(); ++k)
    {
        const double dot = k_fractional[k] * lattice_R;
        h += (k_weights[k] * std::exp(std::complex<double>(0, -2 * PI * dot))) *
             (parent ? states[k].parent_h : states[k].output_h);
    }
    return h;
}

CrpaModelSnapshot snapshot_crpa_model(CrpaPreparedModel& model, const MeanField& meanfield)
{
    CrpaModelSnapshot snapshot;
    snapshot.n_spins = meanfield.get_n_spins();
    snapshot.n_spinor = meanfield.get_n_spinor();
    snapshot.n_aos = meanfield.get_n_aos();
    snapshot.n_bands = meanfield.get_n_bands();
    snapshot.fermi = meanfield.get_efermi();
    snapshot.parent_sites = model.windows.parent_sites;
    snapshot.output_sites = model.windows.sites;
    snapshot.output_columns = model.windows.output_columns;
    const int nk = meanfield.get_n_kpoints();
    snapshot.spin_k.resize(snapshot.n_spins);
    for (int spin = 0; spin < snapshot.n_spins; ++spin)
        for (int k = 0; k < nk; ++k)
        {
            CrpaModelK state;
            state.energies.resize(snapshot.n_bands);
            state.occupations.resize(snapshot.n_bands);
            ComplexMatrix wfc(snapshot.n_bands, snapshot.n_spinor * snapshot.n_aos);
            for (int band = 0; band < snapshot.n_bands; ++band)
            {
                state.energies[band] = meanfield.get_eigenvals()[spin](k, band);
                state.occupations[band] = meanfield.get_weight()[spin](k, band) * nk;
                for (int component = 0; component < snapshot.n_spinor; ++component)
                {
                    const auto* block = meanfield.find_wfc(spin, component, k);
                    if (!block)
                        throw LIBRPA_INVALID_ARGUMENT("model snapshot needs all wavefunctions");
                    for (int ao = 0; ao < snapshot.n_aos; ++ao)
                        wfc(band, component * snapshot.n_aos + ao) = (*block)(band, ao);
                }
            }
            const auto& parent_phi = model.windows.parent_orbitals_spin_k[spin][k];
            state.parent_phi = parent_phi;
            state.parent_capture_min = model.windows.parent_capture_min[spin][k];
            // Use the final cached frame which made T(R) and the U density vertices.
            const auto output_phi = model.frames[spin]->build_phi(k);
            const auto& overlap = model.frames[spin]->get_S_k(k);
            if (model.retain_raw_model && spin == 0)
            {
                ComplexMatrix spatial(snapshot.n_aos, snapshot.n_aos);
                for (int i = 0; i < snapshot.n_aos; ++i)
                    for (int j = 0; j < snapshot.n_aos; ++j) spatial(i, j) = overlap(i, j);
                snapshot.overlap.push_back(std::move(spatial));
            }
            state.parent_projector = transpose(conj(wfc) * (overlap * parent_phi), true);
            state.output_projector = transpose(conj(wfc) * (overlap * output_phi), true);
            if (model.retain_raw_model)
            {
                state.coefficients = transpose(wfc, false);
                state.output_phi = output_phi;
            }
            state.parent_h = project_crpa_one_body(state.parent_projector, state.energies);
            state.output_h = project_crpa_one_body(state.output_projector, state.energies);
            state.parent_density = project_crpa_one_body(state.parent_projector, state.occupations);
            state.output_density = project_crpa_one_body(state.output_projector, state.occupations);
            snapshot.spin_k[spin].push_back(std::move(state));
        }
    return snapshot;
}
}  // namespace librpa_int
