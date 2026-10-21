#include "task_helper.h"

librpa_int::MeanField meanfield_for_output(const librpa::MeanfieldData &snapshot)
{
    librpa_int::MeanField result(snapshot.n_spins, snapshot.n_kpoints, snapshot.n_states,
                                 snapshot.n_aos, snapshot.n_spinor);
    result.get_efermi() = snapshot.fermi_ha;
    for (int spin = 0; spin < snapshot.n_spins; ++spin)
        for (int k = 0; k < snapshot.n_kpoints; ++k)
            for (int state = 0; state < snapshot.n_states; ++state)
            {
                const auto i =
                    (static_cast<std::size_t>(spin) * snapshot.n_kpoints + k) * snapshot.n_states +
                    state;
                result.get_eigenvals()[spin](k, state) = snapshot.energies_ha.at(i);
                result.get_weight()[spin](k, state) = snapshot.weighted_occupations.at(i);
            }
    return result;
}
