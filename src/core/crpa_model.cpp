#include "crpa_model.h"

#include <cmath>

#include "../utils/error.h"

namespace librpa_int
{
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

ComplexMatrix rotate_crpa_tensor(const ComplexMatrix& tensor, const ComplexMatrix& left_rotation,
                                 const ComplexMatrix& right_rotation)
{
    const auto pair_transform = [](const ComplexMatrix& rotation)
    {
        if (rotation.nr <= 0 || rotation.nc <= 0 || rotation.nc > rotation.nr)
            throw LIBRPA_INVALID_ARGUMENT("orbital rotation needs a nonempty isometry");
        const auto gram = transpose(rotation, true) * rotation;
        for (int i = 0; i < gram.nr; ++i)
            for (int j = 0; j < gram.nc; ++j)
                if (!std::isfinite(std::abs(gram(i, j))) ||
                    std::abs(gram(i, j) - (i == j ? 1.0 : 0.0)) > 1e-10)
                    throw LIBRPA_INVALID_ARGUMENT("orbital rotation is not an isometry");
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
