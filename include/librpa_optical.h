#pragma once

#include "librpa_handler.h"
#include "librpa_options.h"

#ifdef __cplusplus
extern "C"
{
#endif

    /** One dense real-space block, with interleaved complex values.
     * translation is an integer lattice translation; row_major is zero or one.
     * The caller owns values_ri and retains it until the collective call returns.
     */
    typedef struct
    {
        int translation[3];
        const double *values_ri;
        int row_major;
    } LibrpaLocalizedBlock;

    typedef struct
    {
        int dimension;
        int n_blocks;
        const LibrpaLocalizedBlock *blocks;
    } LibrpaLocalizedOperator;

    /** Original-H/S/r velocity in the loaded complete SCF eigenframe.
     * Operators are supplied on handler rank zero only. Other ranks pass NULL.
     * H is in Hartree, S dimensionless, and position/lattice lengths in Bohr.
     * AO-spinor rows are AO-major, component-minor. Fourier phase is +i k.R.
     * S and its derivative use Hermitian upper-triangle completion; position
     * completion is r_ij=conj(raw_r_ji+i*dS_ji) for i>=j. These are numerical
     * conventions, independent of a producer filename or input-file format.
     * No SCF energies, occupations or coefficients are changed.
     */
    typedef struct
    {
        int n_spins;
        const LibrpaLocalizedOperator *hamiltonian;
        LibrpaLocalizedOperator overlap;
        LibrpaLocalizedOperator position[3];
        /** Required positive values only for strict-2D optical preparation. */
        double inplane_area_bohr2;
        double auxiliary_monopole_norm_squared;
    } LibrpaNativeVelocityInput;

    /** Collective velocity computation and setter; optional analytic head/wing
     * preparation reuses the ordinary numerical API. Options use API values 3/4,
     * never driver preset 6. Optical preparation requires replace_w_head enabled.
     * Input is borrowed only for this call.
     */
    void librpa_prepare_native_velocity(LibrpaHandler *h, const LibrpaOptions *options,
                                        const LibrpaNativeVelocityInput *root_input,
                                        int prepare_optical);

    /** One rank-local wavefunction block, with band rows and spatial AO columns.
     * Values are row-major interleaved complex numbers. For k-BLACS input the
     * rows and columns are the local dimensions of the handler's descriptor.
     */
    typedef struct
    {
        int spin, spinor, kpoint;
        int rows, columns;
        const double *values_ri;
    } LibrpaOpticalWavefunctionBlock;

    /** Borrowed numerical mean field for an optical calculation. Energies and
     * Fermi level are in Hartree. Spectral arrays have [spin][k][band] layout;
     * occupation_weights use the MeanField convention, occupation/n_kpoints.
     * They may be NULL together only for the optional full-source WFC cache.
     * Wavefunction blocks remain rank-local and use zero-based indices.
     */
    typedef struct
    {
        int n_spins, n_kpoints, n_states, n_aos, n_spinor;
        double fermi_energy_ha;
        const double *eigenvalues_ha;
        const double *occupation_weights;
        int n_wavefunction_blocks;
        const LibrpaOpticalWavefunctionBlock *wavefunction_blocks;
    } LibrpaOpticalMeanFieldInput;

    /** Supplied velocity and optional optical eigenframe. All pointers are
     * borrowed until the collective call returns and are supplied on each rank.
     * Velocity layout is [spin][k][axis][band][band], row-major interleaved
     * complex values in Hartree*Bohr. A NULL meanfield uses the handler SCF data.
     * normalize_ibz_occupations selects star-weighted occupations from the
     * handler; set it to zero for already unweighted optical occupations/nk.
     * Optional full-source inputs are paired WFC/velocity data with fractional
     * coordinates [k][3]; the core derives their k-star member mapping.
     * A zero frequency count selects the ordinary numerical frequency grid.
     */
    typedef struct
    {
        int n_spins, n_kpoints, n_states;
        const double *velocity_ri;
        const LibrpaOpticalMeanFieldInput *meanfield;
        int coefficients_are_kblacs_2d;
        int normalize_ibz_occupations;
        const LibrpaOpticalMeanFieldInput *full_source_meanfield;
        const double *full_source_velocity_ri;
        const double *full_source_kfrac;
        int n_frequencies;
        const double *frequencies_ha;
        double inplane_area_bohr2;
        double auxiliary_monopole_norm_squared;
    } LibrpaSuppliedVelocityInput;

    /** Collective analytic head/wing preparation from numerical velocity data.
     * Options use numerical values 3/4/53 and require replace_w_head enabled.
     * The loaded SCF eigenvalues, occupations and wavefunctions are retained.
     */
    void librpa_prepare_supplied_velocity(LibrpaHandler *h, const LibrpaOptions *options,
                                          const LibrpaSuppliedVelocityInput *input, int need_wing);

    typedef struct LibrpaOpticalResult LibrpaOpticalResult;

    /** Rank-local owning snapshot for retrieving numerical optical outputs.
     * Complex arrays are row-major, interleaved real/imaginary. Velocity layout
     * is [spin][k][axis][row][column]; head/wing are [frequency][row][column].
     * Pointers remain valid until librpa_delete_optical_result is called.
     */
    typedef struct
    {
        int n_spins, n_kpoints, n_states, n_frequencies;
        int head_rows, head_columns, wing_rows, wing_columns;
        const double *velocity_ri;
        const double *frequencies_ha;
        const double *head_ri;
        const double *wing_ri;
    } LibrpaOpticalView;

    LibrpaOpticalResult *librpa_get_optical_result(const LibrpaHandler *h);
    LibrpaOpticalView librpa_get_optical_view(const LibrpaOpticalResult *result);
    void librpa_delete_optical_result(LibrpaOpticalResult *result);

#ifdef __cplusplus
}
#endif
