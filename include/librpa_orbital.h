#pragma once
#include "librpa_handler.h"

#ifdef __cplusplus
extern "C"
{
#endif
    typedef struct LibrpaOrbitalResult LibrpaOrbitalResult;

    typedef struct
    {
        int n_spins, n_spinor, n_kpoints, n_bands, n_aos, n_groups;
        double fermi_ha, highest_occupied_ha;
        double min_overlap_eigenvalue, max_norm_error, max_sum_error, max_sqrt_residual;
    } LibrpaOrbitalInfo;

    typedef struct
    {
        int atom_index, type_index, l;
        const char* channel;  // s, p, eg, t2g, or lN; global Cartesian axes.
    } LibrpaOrbitalGroup;

    typedef struct
    {
        double energy_ha, occupation, kweight, norm;
        const double* mulliken_ao;
        const double* lowdin_ao;
        const double* mulliken_group;
        const double* lowdin_group;
    } LibrpaOrbitalState;

    /** Collective orbital analysis from loaded eigenstates and spatial S(k).
     * S is packed [k][AO row][AO column][real/imag], replicated on every rank.
     * Requires complete scalar or two-component spinor coefficients, shell metadata, explicit basis
     * conventions, and normalized loaded k weights. Evaluates S^(1/2), never H.
     * Mulliken = Re[c* (Sc)]; Lowdin = |S^(1/2)c|^2/(c-dagger S c).
     * For spinors, sum both components in each spatial AO and normalize by their
     * combined S norm. n_spins is one; no extra spin-degeneracy factor is added.
     * All radial shells contribute to atom/l groups, with d split into eg/t2g.
     * Does not alter the handler, eigenstates, cRPA windows or response weights.
     * Each rank owns its result and must delete it. Views remain valid until deletion.
     */
    LibrpaOrbitalResult* librpa_analyze_orbitals(const LibrpaHandler* h, int n_kpoints, int n_aos,
                                                 const double* overlap_k_ri, double residual_tol);
    void librpa_delete_orbital_result(LibrpaOrbitalResult* result);
    void librpa_orbital_info(const LibrpaOrbitalResult* result, LibrpaOrbitalInfo* info);
    void librpa_orbital_group(const LibrpaOrbitalResult* result, int group,
                              LibrpaOrbitalGroup* view);
    void librpa_orbital_state(const LibrpaOrbitalResult* result, int spin, int k, int band,
                              LibrpaOrbitalState* view);

    /** Local Gaussian PDOS from an existing result; no collective or files.
     * Grid is absolute Hartree and sigma is Hartree. Output arrays, supplied by
     * the caller, each have n_spins*n_energies*n_groups doubles in [spin][E][group]
     * order and units states/(Ha cell stored-spin-channel). Uses loaded k weights;
     * applies neither occupation nor an extra spin-degeneracy factor.
     */
    void librpa_orbital_pdos(const LibrpaOrbitalResult* result, int n_energies,
                             const double* energy_grid_ha, double sigma_ha, double* mulliken,
                             double* lowdin);
#ifdef __cplusplus
}
#endif
