#pragma once
#include <stddef.h>

#include "librpa_handler.h"
#include "librpa_optical.h"
#include "librpa_options.h"
#ifdef __cplusplus
extern "C"
{
#endif

    typedef struct LibrpaQsgwTask LibrpaQsgwTask;
    /** All energies are Hartree unless the field explicitly says eV. Strings and
     * numeric inputs are copied by the call receiving them. Checkpoints preserve
     * the immutable input eigenframe; original producer inputs remain required.
     * Keep the handler alive until the task is deleted. After preparation keep
     * its SCF eigenframe/grid fixed; finish band/Coulomb/velocity input before
     * run. Do not modify or reenter the handler during run or its callbacks. Call prepare and set_input collectively before the single run call. */
    typedef struct
    {
        int min_iterations, max_iterations, unoccupied_keep, cut_mode;
        double mixing_beta, cut_shift_ha, convergence_tolerance_ev;
        const char *mixer, *velocity_mode, *restart_file, *checkpoint_file;
        int output_matrices, output_qp, export_hamiltonian;
        int state_low, state_high, interband_velocity;
    } LibrpaQsgwControls;
    typedef struct
    {
        int spin, kpoint, dimension, row_major;
        const double* values_ri;
    } LibrpaQsgwMatrix;
    typedef struct
    {
        int n_vxc, n_vxc_band;
        const LibrpaQsgwMatrix *vxc, *vxc_band;
        int vxc_in_ao_basis;
        const LibrpaNativeVelocityInput* native_operators;
        double inplane_area_bohr2, auxiliary_monopole_norm_squared;
    } LibrpaQsgwInput;
    typedef enum
    {
        LIBRPA_QSGW_ITERATIONS,
        LIBRPA_QSGW_EIGENVALUES,
        LIBRPA_QSGW_MATRICES,
        LIBRPA_QSGW_QUASIPARTICLES,
        LIBRPA_QSGW_ENERGY_QP,
        LIBRPA_QSGW_GW_BAND,
        LIBRPA_QSGW_KS_BAND,
        LIBRPA_QSGW_EXX_BAND,
        LIBRPA_QSGW_QSGW_BAND,
        LIBRPA_QSGW_HAMILTONIAN,
        LIBRPA_QSGW_HEAD_HAMILTONIAN
    } LibrpaQsgwOutputKind;
    /** Borrowed output on rank zero only, valid during the callback. spin is
     * zero-based; text_size excludes the terminating zero. The callback presents
     * output only and must not mutate or recursively use the handler/task. */
    typedef struct
    {
        LibrpaQsgwOutputKind kind;
        int iteration, spin, append;
        const char* text;
        size_t text_size;
        const LibrpaLocalizedOperator* hamiltonian;
    } LibrpaQsgwOutput;
    typedef int (*LibrpaQsgwOutputSink)(const LibrpaQsgwOutput*, void* user_data);
    /** Defaults: min/max=1/10, unoccupied_keep=10, uncut, mixer=none,
     * beta=0.2 (used only for single-k linear mixing), shift=2 eV,
     * convergence=1e-4 eV, reference velocity, QP/HR outputs on, matrices off,
     * restart=none and checkpoint=qsgw_checkpoint.bin. Indices are zero-based. */
    void librpa_init_qsgw_controls(LibrpaQsgwControls*);
    LibrpaQsgwTask* librpa_create_qsgw_task(LibrpaHandler*, const LibrpaOptions*,
                                            const LibrpaQsgwControls*);
    /** Collective validation before ordinary numerical input is installed. */
    void librpa_prepare_qsgw_task(LibrpaQsgwTask*);
    /** Collective. Vxc and native operators are supplied on root; ignored elsewhere. */
    void librpa_set_qsgw_input(LibrpaQsgwTask*, const LibrpaQsgwInput* root_input);
    void librpa_set_qsgw_output_sink(LibrpaQsgwTask*, LibrpaQsgwOutputSink, void* user_data);
    void librpa_run_qsgw_task(LibrpaQsgwTask*);
    void librpa_delete_qsgw_task(LibrpaQsgwTask*);
#ifdef __cplusplus
}
#endif
