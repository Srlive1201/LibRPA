#pragma once

#include <stddef.h>

#include "librpa_handler.h"

#ifdef __cplusplus
extern "C"
{
#endif

    /** Borrowed GW energies for presentation. All energies use Hartree.
     * vxc_ha is [spin][loaded k][all states]; vexx_ha and sigc_ri are
     * [spin][loaded k][calculated state], with interleaved real/imaginary sigc.
     * Output k coordinates are [output k][xyz]. A null output_to_input_kpoint
     * selects the loaded k points in order. The state interval is half open.
     */
    typedef struct
    {
        int band_path;
        int n_spins, n_kpoints, n_states;
        int i_state_low, i_state_high;
        int n_output_kpoints;
        const double *fractional_k;
        const int *output_to_input_kpoint;
        const double *vxc_ha, *vexx_ha, *sigc_ri;
        int n_kpoint_weights;
        const double *kpoint_weights;
    } LibrpaGwOutputInput;

    /** Text is borrowed for this call only. Return zero on success. */
    typedef int (*LibrpaGwTextSink)(void *userdata, const char *text, size_t length);

    /** Noncollective presentation using the loaded mean-field state. Call on the
     * output rank only, after collecting the supplied energy arrays there. Inputs
     * remain valid until return; no handler state or occupations are changed.
     * Band analysis uses the ordinary LibRPA logger. Table calls invoke the sink
     * once with the existing energy_qp or per-spin GW band file format. A sink
     * must not mutate or reenter the handler. These C/C++ helpers have no Fortran
     * wrapper.
     */
    void librpa_print_gw_band_analysis(LibrpaHandler *h, const LibrpaGwOutputInput *input);
    void librpa_write_gw_energy_table(LibrpaHandler *h, const LibrpaGwOutputInput *input,
                                      LibrpaGwTextSink sink, void *userdata);
    void librpa_write_gw_band_spin_table(LibrpaHandler *h, const LibrpaGwOutputInput *input,
                                         int spin, LibrpaGwTextSink sink, void *userdata);

#ifdef __cplusplus
}
#endif
