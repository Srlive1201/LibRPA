#pragma once

#include "librpa_handler.h"
#include "librpa_options.h"
#include "librpa_orbital.h"
#include "librpa_qsgw.h"

#ifdef __cplusplus
extern "C"
{
#endif

    /** Numerical cRPA input. Strings and arrays remain owned by the caller.
     * Initialize with librpa_init_crpa_input before assigning fields. All input
     * values and array contents must be replicated identically on every rank of
     * the handler communicator and remain valid until the collective call returns.
     * Energy windows use inclusive absolute-Hartree lower/upper pairs.
     * B selects the KS space used to project atomic trials. Provide exactly one
     * B window array or explicit band list. With use_orbital_weights OFF, also
     * provide an independent A selector for both branches of Pd. With it ON,
     * provide no A selector: the output frame defines endpoint probabilities.
     * Species labels follow the zero-based type indices passed to librpa_set_atoms.
     * Scalar output uses each spin's parent frame; SOC uses one full spinor frame.
     */
    typedef struct
    {
        int n_species;
        const char* const* species_labels;
        const char* correlated_species;
        const char* ligand_species;
        const char* parent_orbitals;
        const char* output_orbitals;
        int n_response_edges;
        const double* response_windows_ha;
        int n_orbital_edges;
        const double* orbital_windows_ha;
        /** Zero-based indices, shared by all spin/k points. Each list replaces
         * its respective energy window; set the unused count to zero. */
        int n_response_bands;
        const int* response_bands;
        int n_orbital_bands;
        const int* orbital_bands;
        int n_kpoints;
        int n_aos;
        /** Spatial overlap [k][AO row][AO column][real/imag], full k grid.
         * SOC lifts this internally to I_spinor tensor S; n_aos stays spatial. */
        const double* overlap_k_ri;
        double gram_abs;
        double gram_rel;
        double gram_cond_max;
        double s_abs;
        double s_rel;
        double s_cond_max;
        double residual_tol;
        double t_roundtrip_tol;
        /** Off by default: retain C/S/Phi in the owned result for arbitrary AO reprojection.
         * P/H/density and the state energies/occupations are always retained. */
        LibrpaSwitch retain_raw_model;
        /** OFF: binary A windows/lists. ON: probabilities from the actual output
         * frame, w_p=sum_a |Phi_a^dagger S C_p|^2/(C_p^dagger S C_p).
         * Pd weights each original transition by w_p*w_q. No A window/list
         * may be supplied; B and the output orbitals define the projector. */
        LibrpaSwitch use_orbital_weights;
        /** Optional previous metric-orthonormal parent Phi, packed as
         * [spin][k][component-major AO][parent column][real/imag]. A positive
         * n_previous_parent enables transport through the declared current B.
         * The caller owns the reference; the original KS states are unchanged.
         * A zero count retains the original atomic-trial construction. */
        int n_previous_parent;
        const double* previous_parent_phi_ri;
        /** Minimum retained probability of every previous parent direction.
         * Default 0.99. This is an additional subspace-retention gate, not a
         * replacement or relaxation of the Gram/metric rank checks. */
        double tracking_min_capture;
    } LibrpaCrpaInput;

    typedef struct LibrpaCrpaResult LibrpaCrpaResult;

    /** Owned snapshot dimensions and static geometry; views survive handler changes.
     * All indices are zero based. Geometry uses Bohr and fractional k coordinates.
     * n_parent/n_output include both trial components for SOC.
     */
    typedef struct
    {
        int n_spins, n_spinor, n_kpoints, n_bands, n_aos, n_atoms;
        int n_parent, n_output, n_parent_sites, n_output_sites;
        double fermi_ha;
        const double* lattice_bohr;   /* [3][3], row lattice vectors */
        const double* positions_bohr; /* [atom][3] */
        const int* atom_types;
        const int* ao_counts;       /* spatial AO count per atom */
        const double* k_fractional; /* [k][3] */
        const double* k_weights;    /* sum to one */
        const int* output_columns;  /* indices in joint parent frame */
        int n_shells;
        const int* shell_offsets;    /* [atom+1], indexing shell_l */
        const int* shell_l;          /* radial shells in producer AO order */
        const int* basis_convention; /* phase, ratom, order, coeff_negative, coeff_positive */
    } LibrpaCrpaModelInfo;

    /** Row-major complex views packed as [real,imag]; energies use Hartree.
     * C is [n_spinor*n_aos][n_bands], S is spatial [n_aos][n_aos].
     * Phi is [n_spinor*n_aos][n_orbitals], P is [n_orbitals][n_bands].
     * H is [n_orbitals][n_orbitals]; parent and output have separate dimensions.
     * Occupations are physical per-state values, without the k-point weight.
     * C/S/output-Phi pointers are null unless retain_raw_model was enabled.
     * Parent Phi is always retained for numerical subspace transport.
     */
    typedef struct
    {
        const double* energies_ha;
        const double* occupations;
        const double* coefficients_ri;
        const double* overlap_ri;
        const double* parent_phi_ri;
        const double* output_phi_ri;
        const double* parent_projector_ri;
        const double* output_projector_ri;
        const double* parent_h_ha_ri;
        const double* output_h_ha_ri;
        const double* parent_density_ri;
        const double* output_density_ri;
        double parent_capture_min;
    } LibrpaCrpaModelK;

    void librpa_crpa_model_info(const LibrpaCrpaResult* result, LibrpaCrpaModelInfo* info);
    void librpa_crpa_model_k(const LibrpaCrpaResult* result, int spin, int kpoint,
                             LibrpaCrpaModelK* view);
    void librpa_crpa_model_site(const LibrpaCrpaResult* result, int parent, int site, int* atom,
                                int* column_start, int* n_orbitals);
    const char* librpa_crpa_model_label(const LibrpaCrpaResult* result, int parent, int site,
                                        int orbital);
    /** Weighted onsite H and density in one site's actual frame, supplied output buffers [M][M][2].
     */
    void librpa_crpa_model_local(const LibrpaCrpaResult* result, int parent, int spin, int site,
                                 double* h_ha_ri, double* density_ri);
    /** Full parent/output hopping at integer lattice R, [M][M][2], Hartree.
     * Uses sum_k weight_k exp(-i 2*pi k.R) H(k). R need not be a stored BvK representative. */
    void librpa_crpa_model_hopping(const LibrpaCrpaResult* result, int parent, int spin,
                                   const int* lattice_R, double* h_ha_ri);
    /** Noncollective numerical helpers; arrays are row-major complex pairs.
     * P=[orbital][band], diagonal uses the caller's energy or occupation convention. */
    void librpa_crpa_project_one_body(int n_orbitals, int n_bands, const double* projector_ri,
                                      const double* diagonal, double* operator_ri);
    /** Construct P from an explicit metric-orthonormal frame without a KS solve.
     * C=[component-major AO][band], Phi=[component-major AO][orbital], S is spatial.
     * Returns P=[orbital][band]. The orbital construction/window is chosen by the caller. */
    void librpa_crpa_project_orbitals(int n_aos, int n_spinor, int n_bands, int n_orbitals,
                                      const double* coefficients_ri, const double* overlap_ri,
                                      const double* phi_ri, double* projector_ri);
    /** P_new=R^dagger P_old for Phi_new=Phi_old R; use the same R on the corresponding U legs. */
    void librpa_crpa_rotate_projector(int n_bands, int n_old, int n_new, const double* projector_ri,
                                      const double* rotation_ri, double* rotated_projector_ri);
    /** Phi_new=Phi_old R, R=[old][new] is an isometry; rotate both density legs.
     * No screening response is changed. Use the matching rotated projector/H. */
    void librpa_crpa_rotate_tensor(int n_old, int n_new, const double* tensor_ri,
                                   const double* left_rotation_ri, const double* right_rotation_ri,
                                   double* rotated_ri);
    /** Density U_ab,cd -> fermionic V_acbd: V_ijkl=U_ik,jl.
     * V multiplies (1/2) c_i^dagger c_j^dagger c_l c_k, with the stored spin-frame pair.
     * No antisymmetrization or Kanamori/Slater approximation is performed. */
    void librpa_crpa_tensor_to_fermion(int n_orbitals, const double* density_ri,
                                       double* fermion_ri);

    /** Ordered-column averages of U_mm,mm, U_mm,nn and U_mn,nm (m!=n).
     * Input is [M][M][M][M][2]; output is three complex pairs in the same energy unit.
     * Off-diagonal means are zero for M=1. Spinor columns use the same ordered means,
     * without a scalar-shell Hund interpretation or an added degeneracy factor. */
    void librpa_crpa_tensor_averages(int n_orbitals, const double* tensor_ri,
                                     double* averages_ha_ri);

    typedef struct LibrpaCrpaContinuedResult LibrpaCrpaContinuedResult;
    typedef struct
    {
        int atom_index, spin_left, spin_right, n_orbitals;
        const double* values_ha_ri;
        const int* status;              /* [M^4], zero means numerical diagnostics pass */
        const double* subset_spread_ha; /* [M^4] */
        double max_source_residual_ha, max_holdout_residual_ha;
    } LibrpaCrpaContinuedTensor;
    /** Native Padé of U-v using original z=i*nu; restores v, no chemical-potential shift.
     * Queries are [real,imag] in Hartree: z=0 for static extrapolation, or Im(z)>0.
     * n_parameters=-1 uses all nodes. tolerance_ha controls numerical diagnostics.
     * Noncollective; returns independently owned values, not an added response node.
     * Status bits: 1 primary failure, 2 subset failure, 4 excessive subset spread,
     * 8 excessive source/holdout residual, 16 fewer than four nodes, 32 nonfinite value.
     * No passivity, causality, basis or mesh convergence certificate is implied. */
    LibrpaCrpaContinuedResult* librpa_continue_crpa(const LibrpaCrpaResult* result, int n_queries,
                                                    const double* queries_ha_ri, int n_parameters,
                                                    double tolerance_ha);
    void librpa_delete_crpa_continued_result(LibrpaCrpaContinuedResult* result);
    int librpa_crpa_continued_size(const LibrpaCrpaContinuedResult* result);
    void librpa_crpa_continued_tensor(const LibrpaCrpaContinuedResult* result, int tensor,
                                      int query, LibrpaCrpaContinuedTensor* view);

    /** Read-only tensor view. Arrays remain valid until the result is deleted.
     * Arrays contain M^4 complex entries, each packed as real/imag doubles.
     * Matrix row is the ordered pair (a,b), column is (c,d): a*M+b, c*M+d.
     * Energies and interactions are Hartree; frequency is positive imaginary nu.
     */
    typedef struct
    {
        int atom_index;
        int spin_left;
        int spin_right;
        int n_orbitals;
        double frequency_ha;
        const double* bare_ri;
        const double* partially_screened_ri;
        /** Null for U-only results; no fully screened interaction was computed. */
        const double* fully_screened_ri;
    } LibrpaCrpaTensor;

    void librpa_init_crpa_input(LibrpaCrpaInput* input);

    /** Read full fractional k coordinates in the current KS eigenvector order. */
    void librpa_get_crpa_kgrid(const LibrpaHandler* h, int n_kpoints, double* fractional_k);
    /** Prepare current local frames, weights and one-body model only. Collective.
     * Builds no response, W or self-energy; returns zero interaction tensors.
     * Useful for cheap parent transport between selected screening evaluations. */
    LibrpaCrpaResult* librpa_prepare_crpa_model(const LibrpaHandler* h,
                                                const LibrpaOptions* options,
                                                const LibrpaCrpaInput* input);

    /** Collective calculation on the handler communicator; no input/output files.
     * All ranks must call with identical options and replicated LibrpaCrpaInput.
     * Supports scalar or two-component spinor full-grid LIBRI with replicated KS eigenvectors and
     * native positive minimax nodes. All EXX/GW/RPA symmetry switches must be disabled.
     * No static node or analytic continuation is added. The result is independently owned on each
     * rank; delete it on every rank.
     */
    LibrpaCrpaResult* librpa_compute_crpa_window(const LibrpaHandler* h,
                                                 const LibrpaOptions* options,
                                                 const LibrpaCrpaInput* input);

    /** Build the ordinary GW self-energy and cRPA U for the same entering state.
     * P0 is built once by GW and borrowed by cRPA; only Pd and constrained
     * screening are added. GW P0, optical data, W and self-energy conventions
     * remain unchanged. Requires full-grid LIBRI without symmetry or a SigC
     * file-read shortcut. Input has the same contract as the standalone call.
     * The owned result contains bare/U; fully_screened_ri is null.
     * cRPA screens and projects with full Coulomb independently of the GW
     * eps/Wc/EXX Coulomb flags, which keep their original GW meaning.
     */
    LibrpaCrpaResult* librpa_build_g0w0_sigma_crpa(LibrpaHandler* h, const LibrpaOptions* options,
                                                   const LibrpaCrpaInput* input);

    /** Response-only bare/U calculation for the handler's current E/C/f.
     * Builds P0 and Pd and screens Pr, without full-W screening or self-energy.
     * Suitable for the final updated QSGW state. Uses current velocity data
     * when optical correction is enabled. Returns an independently owned
     * result whose fully_screened_ri is null. The cRPA solve uses full Coulomb
     * regardless of the independent GW Coulomb flags.
     */
    LibrpaCrpaResult* librpa_compute_crpa_u_window(const LibrpaHandler* h,
                                                   const LibrpaOptions* options,
                                                   const LibrpaCrpaInput* input);

    /** Noncollective save of this result's parent frame in the native tracking
     * checkpoint format. The caller selects the writing rank. Energy context
     * may be empty for a state export, or include every spin/k/band pole. */
    void librpa_write_crpa_tracking_checkpoint(const LibrpaCrpaResult* result, const char* filename,
                                               int iteration, const char* definition,
                                               int n_energies, const double* energies_ha);

    typedef struct
    {
        int interval, track_orbitals, output_analysis;
        const char* tracking_definition;
    } LibrpaQsgwCrpaControls;

    /** Borrowed numerical output on every rank, valid only during the call.
     * The sink may query results and write files, but must not delete these
     * views or recursively use the handler/task. Return zero on success. */
    typedef int (*LibrpaCrpaOutputSink)(const LibrpaCrpaResult*, const LibrpaOrbitalResult*,
                                        int state_index, int gw_update, int final_state,
                                        void* user_data);
    void librpa_init_qsgw_crpa_controls(LibrpaQsgwCrpaControls* controls);
    /** Attach after librpa_prepare_qsgw_task, before run. Input and controls
     * are copied; the sink and user_data must remain valid until run returns. */
    void librpa_set_qsgw_crpa(LibrpaQsgwTask* task, const LibrpaCrpaInput* input,
                              const LibrpaQsgwCrpaControls* controls, LibrpaCrpaOutputSink sink,
                              void* user_data);

    void librpa_delete_crpa_result(LibrpaCrpaResult* result);
    int librpa_crpa_result_size(const LibrpaCrpaResult* result);
    void librpa_crpa_result_tensor(const LibrpaCrpaResult* result, int index,
                                   LibrpaCrpaTensor* tensor);

    /** Retrieve A (response_window != 0) or B selected original-KS band indices.
     * Band indices are zero based and the returned array is owned by result.
     */
    void librpa_crpa_result_bands(const LibrpaCrpaResult* result, int response_window, int spin,
                                  int kpoint, int* count, const int** bands);
    /** Actual endpoint weights [band], owned by result. Binary or fractional,
     * depending on the input mode. Weighted mode has no binary A band list. */
    void librpa_crpa_result_weights(const LibrpaCrpaResult* result, int spin, int kpoint,
                                    int* count, const double** weights);

    /** Trial label for a local output column; owned by result until deletion.
     * SOC suffix up/down identifies a trial, not a conserved spin quantum number.
     * SOC tensor spin_left/right are channel zero; spinor indices are a,b,c,d. */
    const char* librpa_crpa_result_orbital_label(const LibrpaCrpaResult* result, int atom,
                                                 int orbital);

#ifdef __cplusplus
}
#endif
