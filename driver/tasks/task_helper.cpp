#include "task_helper.h"

#include <stdexcept>

std::size_t meanfield_output_index(const librpa::MeanfieldData &snapshot, int spin, int k,
                                   int state)
{
    return (static_cast<std::size_t>(spin) * snapshot.n_kpoints + k) * snapshot.n_states + state;
}

GwOutputView::GwOutputView(const librpa::MeanfieldData &snapshot, bool band_path,
                           const std::vector<librpa_int::Vector3_Order<double>> &kfrac,
                           const std::vector<int> &mapping,
                           const std::vector<librpa_int::matrix> &potentials,
                           const std::vector<double> &vexx,
                           const std::vector<std::complex<double>> &sigc, int state_low,
                           int state_high, const std::vector<double> &weights)
{
    const auto count =
        static_cast<std::size_t>(snapshot.n_spins) * snapshot.n_kpoints * (state_high - state_low);
    if (state_low < 0 || state_high <= state_low || state_high > snapshot.n_states ||
        vexx.size() != count || sigc.size() != count ||
        potentials.size() != static_cast<std::size_t>(snapshot.n_spins) ||
        (!mapping.empty() && mapping.size() != kfrac.size()))
        throw std::invalid_argument("Incomplete GW output buffers");
    for (const auto &k : kfrac) fractional_k.insert(fractional_k.end(), {k.x, k.y, k.z});
    for (const auto &potential : potentials)
    {
        if (potential.nr != snapshot.n_kpoints || potential.nc != snapshot.n_states)
            throw std::invalid_argument("Invalid GW output potential shape");
        for (int k = 0; k < snapshot.n_kpoints; ++k)
            for (int state = 0; state < snapshot.n_states; ++state)
                vxc.push_back(potential(k, state));
    }
    input.band_path = band_path;
    input.n_spins = snapshot.n_spins;
    input.n_kpoints = snapshot.n_kpoints;
    input.n_states = snapshot.n_states;
    input.i_state_low = state_low;
    input.i_state_high = state_high;
    input.n_output_kpoints = kfrac.size();
    input.fractional_k = fractional_k.data();
    input.output_to_input_kpoint = mapping.empty() ? nullptr : mapping.data();
    input.vxc_ha = vxc.data();
    input.vexx_ha = vexx.data();
    input.sigc_ri = reinterpret_cast<const double *>(sigc.data());
    input.n_kpoint_weights = weights.size();
    input.kpoint_weights = weights.empty() ? nullptr : weights.data();
}
