#pragma once

#include <string>
#include <vector>

#include "atomic_basis.h"
#include "correlated_subspace.h"
#include "meanfield.h"

namespace librpa_int
{

struct CrpaWindowOptions
{
    std::string correlated_species;
    std::string ligand_species;
    //! Parent frame and output subset: s, p, sp, d, dp, eg, or t2g, in the Cartesian frame.
    //! d|dp and eg|dp select columns after joint d+p polar normalization.
    std::string parent_orbitals = "d";
    std::string output_orbitals = "d";
    bool use_orbital_weights = false;
    //! Inclusive absolute-Hartree bounds, stored as lower/upper pairs.
    std::vector<double> response_windows_ha;
    std::vector<double> orbital_windows_ha;
    //! Explicit zero-based KS indices, identical for all spin/k points.
    //! Each list is mutually exclusive with the corresponding energy windows.
    std::vector<int> response_bands;
    std::vector<int> orbital_bands;
    //! Optional previous parent Phi[spin][k], transported into the declared B.
    //! This tracks the complete parent (including ligands), never just output d.
    std::vector<std::vector<ComplexMatrix>> previous_parent_spin_k;
    double tracking_min_capture = 0.99;
    double gram_abs = 1.0e-12;
    double gram_rel = 1.0e-12;
    double gram_cond_max = 1.0e12;
    double s_abs = 1.0e-10;
    double s_rel = 1.0e-12;
    double s_cond_max = 1.0e10;
    double residual_tol = 1.0e-10;
    double t_roundtrip_tol = 1.0e-10;
};

struct CrpaWindowBuildResult
{
    std::vector<SiteOrbitalGroup> sites;
    //! Joint parent frame before output selection, including ligand sites.
    std::vector<SiteOrbitalGroup> parent_sites;
    std::vector<int> output_columns;
    std::vector<std::vector<ComplexMatrix>> parent_orbitals_spin_k;
    //! Minimum retained probability of a previous parent direction; 1 initially.
    std::vector<std::vector<double>> parent_capture_min;
    //! Phi[channel][k], with component-major spinor AO rows for SOC.
    std::vector<std::vector<ComplexMatrix>> orbitals_spin_k;
    //! Original-KS state indices [spin][k][selected band], independently selected.
    std::vector<std::vector<std::vector<int>>> response_bands;
    std::vector<std::vector<std::vector<int>>> orbital_bands;
    //! Endpoint probabilities from the actual output frame; empty in binary mode.
    std::vector<std::vector<std::vector<double>>> response_weights;
};

std::vector<std::vector<std::vector<int>>> crpa_bands_in_energy_windows(
    const MeanField& meanfield, const std::vector<double>& edges_ha);

//! Lift a spatial overlap to component-major I_spinor tensor S (one or two components).
ComplexMatrix crpa_spinor_overlap(const ComplexMatrix& overlap, int n_spinor);

/*! Build spin-resolved local orbitals without changing KS energies or occupations.
 * Atomic types index species_labels; unrelated species remain spectators.
 * Only the first radial shell for each selected angular momentum is used.
 * Spatial S(k) and the scalar/spinor eigenvectors share the supplied basis convention.
 * SOC uses one jointly normalized spinor frame, not separate physical-spin frames.
 * Response window A never participates in construction of output window B.
 */
CrpaWindowBuildResult build_crpa_window_input(
    const AtomicBasis& basis, const BasisConvention& convention, const std::vector<int>& atom_types,
    const std::vector<std::string>& species_labels, const std::vector<ComplexMatrix>& overlap_k,
    const std::vector<Vector3_Order<double>>& kfrac_list, const MeanField& meanfield,
    const CrpaWindowOptions& options);

}  // namespace librpa_int
