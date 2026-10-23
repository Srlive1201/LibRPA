#pragma once

/*!
 * @file correlated_subspace.h
 * @brief Orthonormal local orbitals and their lattice Fourier transform.
 *
 * Atomic trials are first projected into a supplied KS window B:
 *   W(k) = C_B(k) C_B(k)^dagger S(k) T.
 * The local frame is the metric polar factor:
 *   Phi(k) = W(k) [W(k)^dagger S(k) W(k)]^(-1/2).
 * This class constructs output orbitals; it does not select response bands.
 * Construct a separate instance for each scalar output spin frame. SOC uses
 * one instance with component-major spinor AO rows and I_spinor tensor S.
 *
 * Storage convention (verified against production code):
 *   MeanField::wfc is stored as bands×NAO, where wfc(ib, iao) = C_{iao, ib},
 *   i.e. wfc = C^T (transpose, NOT conjugate-transpose).  This is confirmed
 *   by src/api/input.cpp (raw copy into wfc.c), driver/reader_eigenvec.cpp
 *   (wfc_index with ib major, iw minor), and src/core/meanfield.cpp which
 *   builds the AO Green kernel as transpose(wfc, false) * conj(wfc_scaled).
 *
 * Therefore:
 *   C† = conj(wfc)             (element-wise conjugate of the bands×NAO matrix)
 *   Q  = conj(wfc) * S * Φ     (bands × n_corr)
 *   C†SC = conj(wfc) * S * transpose(wfc, false)
 *
 * Fourier convention (eq. v5-T-fourier):
 *   Φ(k) = Σ_R exp(+i2π k_frac·R) T(R)
 *   T(R) = (1/Nk) Σ_k exp(-i2π k_frac·R) Φ(k)
 *
 * Rank-deficient or ill-conditioned O(k) throws; no silent truncation.
 * Complex S/W/C/Φ/Q supported; spinor rows use the same metric algebra.
 */

#include <complex>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "../math/complexmatrix.h"
#include "../math/vector3_order.h"
#include "../utils/constants.h"
#include "../utils/error.h"

namespace librpa_int
{

/*! Project atomic trial columns into an explicit KS band window:
 * W_window = C_window C_window^dagger S W. The stored wfc is C^T.
 * The result is deliberately not normalized; build_phi performs the existing
 * S-metric rank check and polar normalization. No KS energies are changed.
 * This constructs output orbitals; it does not choose a response exclusion.
 */
ComplexMatrix project_trial_orbitals_to_band_window(const ComplexMatrix& wfc,
                                                    const ComplexMatrix& S,
                                                    const ComplexMatrix& trials,
                                                    const std::vector<int>& bands);

/*! Column group for one correlated site within the global trial matrix W. */
struct SiteOrbitalGroup
{
    std::string label;    //!< e.g. "Ni-A", "Ni-B"
    int orb_start;        //!< column offset in W / Φ / Q / T
    int n_orbitals;       //!< e.g. 5 for d shell
    int atom_index = -1;  //!< production atom index; -1 = unset
    std::vector<std::string>
        orbital_labels;  //!< e.g. {"dz²","dxz","dyz","dx²-y²","dxy"}; empty = unchecked
};

/*! Per-k diagnostics for the overlap/Gram part (spin-independent). */
struct SubspaceKptDiagnostics
{
    double gram_lambda_min;    //!< min eigenvalue of O(k)
    double gram_lambda_max;    //!< max eigenvalue of O(k)
    double gram_condition;     //!< λmax / λmin
    double phi_sphi_residual;  //!< ||Φ†SΦ - I||
    double S_lambda_min;       //!< min eigenvalue of S(k)
    double S_lambda_max;       //!< max eigenvalue of S(k)
    double S_condition;        //!< λmax_S / λmin_S
    bool band_basis_complete;  //!< nbands == nao
    int rank;                  //!< numerical rank of O(k)
};

/*! Per-(spin, k) band/Q diagnostics. */
struct SpinKptDiagnostics
{
    double csc_residual;             //!< ||C†SC - I|| for the supplied KS band space
    double q_completeness_residual;  //!< ||Q†Q - I|| for the supplied KS band space
    double qqt_lambda_min;           //!< min eigenvalue of Q†Q (always computed)
    double qqt_lambda_max;           //!< max eigenvalue of Q†Q (always computed)
    double qqt_captured_fraction;    //!< trace(Q†Q) / n_corr, in [0,1] for an S-orthonormal KS band
                                     //!< space
    bool computed;                   //!< whether this (spin,k) has been built
};

/*! Aggregate diagnostics across all k / spin. */
struct SubspaceDiagnostics
{
    double max_gram_condition;
    double min_gram_lambda;
    double max_gram_lambda;
    double max_phi_sphi_residual;
    double max_csc_residual;
    double max_q_completeness_residual;
    double max_t_roundtrip_residual;
    double min_S_lambda;
    double max_S_lambda;
    double max_S_condition;
    int n_kpoints;
    int n_spins;
    int nao;
    int nbands;
    int n_corr;
    int n_sites;
    bool t_computed;
};

/*! Result returned by build_spin_k. */
struct SpinKBuildResult
{
    ComplexMatrix Q;                   //!< bands × n_corr
    SubspaceKptDiagnostics kpt_diag;   //!< spin-independent k-point diagnostics
    SpinKptDiagnostics spin_kpt_diag;  //!< per-(spin,k) band/Q diagnostics
};

class CorrelatedSubspace
{
public:
    /*!
     * @brief Construct with labels, k/R lists and dimensions.
     * @param sites  Site/orbital groups; orb_start must be contiguous and
     *               cover [0, n_corr) without gaps or overlaps.
     * @param kfrac_list  Fractional k-point coordinates (must match Dataset).
     * @param R_list      Lattice translation vectors for T(R).
     * @param nao         Number of atomic orbitals.
     * @param nbands      Number of bands in wfc.
     * @param n_spins     Number of spin channels (1 or 2 for scalar).
     * @param abs_threshold   Absolute eigenvalue floor for O(k).
     * @param rel_threshold   Relative (to λmax) eigenvalue floor for O(k).
     * @param cond_max        Maximum allowed condition number λmax/λmin for O(k).
     * @param s_abs_threshold  Absolute eigenvalue floor for S(k) positive-definiteness.
     * @param s_rel_threshold  Relative (to λmax_S) eigenvalue floor for S(k).
     * @param s_cond_max       Maximum allowed condition number for S(k).
     * @param residual_tol     Residual tolerance for ||Φ†SΦ-I||, ||Q†Q-I||, ||C†SC-I||
     *                         (scaled by sqrt(dimension) internally).
     * @param t_roundtrip_tol  Maximum allowed T(R)↔Φ(k) round-trip Frobenius residual.
     */
    CorrelatedSubspace(const std::vector<SiteOrbitalGroup>& sites,
                       const std::vector<Vector3_Order<double>>& kfrac_list,
                       const std::vector<Vector3_Order<int>>& R_list, int nao, int nbands,
                       int n_spins, double abs_threshold = 1.0e-12, double rel_threshold = 1.0e-12,
                       double cond_max = 1.0e12, double s_abs_threshold = 1.0e-10,
                       double s_rel_threshold = 1.0e-12, double s_cond_max = 1.0e10,
                       double residual_tol = 1.0e-10, double t_roundtrip_tol = 1.0e-10);

    /*! Number of correlated orbitals (total columns of W/Φ). */
    int n_corr() const { return n_corr_; }

    /*! Validate site groups: contiguous, no gaps, no overlaps, no out-of-bounds. */
    static void validate_site_groups(const std::vector<SiteOrbitalGroup>& sites, int n_corr);

    /*! Set S(k) for one k-point. Must be nao×nao Hermitian. Invalidates caches. */
    void set_S_k(int ik, const ComplexMatrix& S);

    /*! Set W(k) for one k-point. Must be nao×n_corr. Invalidates caches. */
    void set_W_k(int ik, const ComplexMatrix& W);

    /*!
     * @brief Validate KS orthogonality and orbital capture for one (spin, k) from bands×NAO wfc.
     *
     * wfc is bands×NAO with wfc(ib, iao) = C_{iao, ib} (i.e. C^T, NOT C†).
     * Computes Q = conj(wfc) * S * Φ for validation only.
     * Throws on rank deficiency or ill-conditioning of O(k).
     * Uses cache-valid fast path: if S/W unchanged since last build_phi(ik),
     * returns the cached Φ without recomputation.
     */
    SpinKBuildResult build_spin_k(int ispin, int ik, const ComplexMatrix& wfc);

    /*!
     * @brief Build Φ(k) = W(k) O^{-1/2}(k) for one k-point.
     * Also computes and stores O(k) and S(k) diagnostics.
     * Performs fail-closed S(k) Hermiticity and positive-definiteness check.
     * Throws on rank deficiency or ill-conditioning.
     * Uses cache-valid fast path: if S/W unchanged, returns cached Φ.
     */
    ComplexMatrix build_phi(int ik);

    /*!
     * @brief Compute real-space T(R) from Φ(k) via Fourier.
     *
     * T(R) = (1/Nk) Σ_k exp(-i2π k_frac·R) Φ(k)
     * Must be called after S and W are set for all k.
     * Fail-closed: verifies R_list/kfrac_list duality and Fourier closure.
     */
    void compute_T_R();

    /*! Get T(R) for one R index. Throws if compute_T_R not called. */
    const ComplexMatrix& get_T(int iR) const;

    /*! Get T(R) columns for one site. Returns nao × n_orbitals_site. */
    ComplexMatrix get_T_site(int iR, int site_index) const;

    /*! Fourier round-trip residual: max ||Φ(k) - Σ_R exp(+ikR) T(R)||. */
    double compute_t_roundtrip_residual() const;

    /*! Access site groups. */
    const std::vector<SiteOrbitalGroup>& get_sites() const { return sites_; }
    //! Read the cached scalar/spinor metric; numerical model export must not duplicate every S(k).
    const ComplexMatrix& get_S_k(int ik) const { return S_k_.at(ik); }

    /*! Access kfrac list. */
    const std::vector<Vector3_Order<double>>& get_kfrac_list() const { return kfrac_list_; }

    /*! Access R list. */
    const std::vector<Vector3_Order<int>>& get_R_list() const { return R_list_; }

    /*! Aggregate diagnostics (valid after all build calls). */
    const SubspaceDiagnostics& get_diagnostics() const { return diag_; }

    /*! Get per-k diagnostics (spin-independent Gram/Φ part). */
    const std::vector<SubspaceKptDiagnostics>& get_kpt_diagnostics() const { return kpt_diag_; }

    /*! Get per-(spin,k) band/Q diagnostics. */
    const std::vector<std::vector<SpinKptDiagnostics>>& get_spin_kpt_diagnostics() const
    {
        return spin_kpt_diag_;
    }

private:
    std::vector<SiteOrbitalGroup> sites_;
    std::vector<Vector3_Order<double>> kfrac_list_;
    std::vector<Vector3_Order<int>> R_list_;

    int nao_;
    int nbands_;
    int n_spins_;
    int n_corr_;
    int n_kpoints_;

    // Gram/O thresholds
    double abs_threshold_;
    double rel_threshold_;
    double cond_max_;

    // S thresholds
    double s_abs_threshold_;
    double s_rel_threshold_;
    double s_cond_max_;

    // Residual thresholds
    double residual_tol_;
    double t_roundtrip_tol_;

    // Per-k data (spin-independent)
    std::vector<ComplexMatrix> S_k_;    //!< nao × nao
    std::vector<ComplexMatrix> W_k_;    //!< nao × n_corr
    std::vector<ComplexMatrix> phi_k_;  //!< nao × n_corr (cached, invalidated on set_S/set_W)
    std::vector<bool> phi_valid_;       //!< whether phi_k_[ik] is up-to-date
    std::vector<SubspaceKptDiagnostics> kpt_diag_;
    std::vector<bool> kpt_diag_valid_;  //!< whether kpt_diag_[ik] is up-to-date

    // Per-(spin,k) data
    std::vector<std::vector<SpinKptDiagnostics>> spin_kpt_diag_;  //!< [ispin][ik]

    // Real-space T(R)
    std::vector<ComplexMatrix> T_R_;  //!< nao × n_corr
    bool T_valid_;                    //!< whether T_R_ is up-to-date

    // Aggregate diagnostics
    SubspaceDiagnostics diag_;

    /*! Invalidate phi, kpt_diag, spin_kpt_diag, and T caches for one k-point.
     *  Immediately recomputes aggregate diagnostics so get_diagnostics()
     *  never returns stale values.
     */
    void invalidate_k(int ik);

    /*! Reset aggregate diagnostics to sentinel values. */
    void reset_diagnostics();

    /*! Diagonalize a Hermitian ComplexMatrix using heev (row-major). */
    void diagonalize_hermitian(ComplexMatrix& mat, std::vector<double>& eigenvalues) const;

    /*! Compute Fourier phase exp(+i 2π k_frac · R) using project TWO_PI. */
    std::complex<double> fourier_phase(const Vector3_Order<double>& kfrac,
                                       const Vector3_Order<int>& R) const;

    /*! Check matrix dimensions and finiteness. */
    void check_matrix(const ComplexMatrix& m, int expected_nr, int expected_nc,
                      const std::string& name) const;

    /*!
     * Check S(k) Hermiticity (full-matrix scale-aware tolerance, diagonal
     * imaginary part) and positive-definiteness (abs + rel thresholds).
     * Records S eigenvalue diagnostics into kpt_diag_[ik].
     * Throws on failure.
     */
    void check_S_k(int ik);

    /*! Compute Frobenius norm of a ComplexMatrix. */
    double frobenius_norm(const ComplexMatrix& m) const;

    /*! Recompute aggregate diagnostics from all valid per-k and per-(spin,k). */
    void recompute_aggregate();

    /*! Verify that kfrac_list and R_list form a dual Fourier pair. Throws if not. */
    void verify_fourier_duality() const;
};

}  // namespace librpa_int
