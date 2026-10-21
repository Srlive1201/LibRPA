/*!
 @file meanfield.h
 @brief Utilities to handle the mean-field starting point for many-body calculation
 */
#ifndef MEANFIELD_H
#define MEANFIELD_H

#include <array>
#include <vector>
#include <map>
#include <utility>
#include "atom.h"
#include "pbc.h"
#include "symmetry_context.h"
#include "symmetry_spin_kernel.h"
#include "../math/matrix.h"
#include "../math/complexmatrix.h"
#include "../math/vector3_order.h"

namespace librpa_int {

//! Object of the meanfield input
/*!
  @note Energies are saved in Hartree unit.
 */
class MeanField
{
private:
    //! number of spin channels
    int n_spins;
    //! number of atomic orbitals
    int n_aos;
    //! number of bands
    int n_states;
    //! number of kpoints
    int n_kpoints;
    //! number of spinor components per wavefunction, e.g. 1 or 2
    int n_spinor;

    // Local dimensions for parallel distribution of eigenvectors
    int n_aos_local;
    int i_ao_start;
    int n_states_local;
    int i_state_start;

    //! eigenvalues, (n_spins, n_kpoints, n_states)
    std::vector<matrix> eskb;
    //! occupation weight, scaled by n_kpoints. (n_spins, n_kpoints, n_states)
    std::vector<matrix> wg;
    //! eigenvector, (n_spins, n_spinor, n_kpoints_local, n_states_local, n_aos_local)
    // TODO: parallelize
    std::map<int, std::map<int, std::map<int, ComplexMatrix>>> wfc;
    //! Fermi energy
    double efermi;
    void resize(int ns, int nk, int nb, int nao, int n_spinor, int st_ib, int nb_local, int st_iao, int nao_local);

public:
    MeanField()
        : n_spins(0),
          n_aos(0),
          n_states(0),
          n_kpoints(0),
          n_spinor(0),
          n_aos_local(0),
          i_ao_start(-1),
          n_states_local(0),
          i_state_start(-1),
          eskb(),
          wg(),
          wfc(),
          efermi(0) {};
    MeanField(int ns, int nk, int nb, int nao, int n_spinor = 1);
    MeanField(int ns, int nk, int nb, int nao, int n_spinor, int st_ib, int nb_local, int st_iao, int nao_local);
    MeanField(int ns, int nk, int nb, int nao, int st_ib, int nb_local, int st_iao, int nao_local)  // backward compability
        : MeanField(ns, nk, nb, nao, 1, st_ib, nb_local, st_iao, nao_local)
    {}
    ~MeanField() {};
    void set(int ns, int nk, int nb, int nao, int n_spinor = 1);
    void set(int ns, int nk, int nb, int nao, int n_spinor, int st_ib, int nb_local, int st_iao, int nao_local);
    void set(int ns, int nk, int nb, int nao, int st_ib, int nb_local, int st_iao, int nao_local)  // backward compability
    {
        set(ns, nk, nb, nao, 1, st_ib, nb_local, st_iao, nao_local);
    }
    MeanField(const MeanField&) = default;
    bool initialized() const { return n_spins > 0 && n_kpoints > 0 && n_states > 0 && n_aos > 0; }
    inline int get_n_bands() const { return n_states; }
    inline int get_n_states() const { return n_states; } // alias
    inline int get_n_spins() const { return n_spins; }
    inline int get_n_kpoints() const { return n_kpoints; }
    inline int get_n_aos() const { return n_aos; }
    inline int get_n_spinor() const { return n_spinor; }
    inline double& get_efermi() { return efermi; }
    inline const double& get_efermi() const { return efermi; }
    std::vector<matrix>& get_eigenvals() { return eskb; }
    const std::vector<matrix>& get_eigenvals() const { return eskb; }
    std::vector<matrix>& get_weight() { return wg; }
    const std::vector<matrix>& get_weight() const { return wg; }
    std::map<int, std::map<int, std::map<int, ComplexMatrix>>>& get_eigenvectors() { return wfc; }
    const std::map<int, std::map<int, std::map<int, ComplexMatrix>>>& get_eigenvectors() const { return wfc; }
    ComplexMatrix* find_wfc(int ispin, int ispinor, int ikpt) noexcept;
    const ComplexMatrix* find_wfc(int ispin, int ispinor, int ikpt) const noexcept;
    double get_E_min_max(double& emin, double& emax) const;
    double get_band_gap() const;
    //! Highest-energy occupied state for a spin. Use ikpt < 0 to search all k-points.
    //! occupation_tol is in unscaled occupation units; stored k-point weights are
    //! compared against occupation_tol / n_kpoints.
    //! Returns {-1, -1} when no occupied state is found.
    std::pair<int, int> find_highest_occupied_state(int ispin, int ikpt = -1,
                                                    double occupation_tol = 1.0e-8) const;
    //! Largest state index whose state and all lower-indexed states are below energy.
    //! Returns -1 if no state satisfies the condition.
    int get_max_state_below_energy(double energy) const;
    //! Smallest state index whose state and all higher-indexed states are above energy.
    //! Returns n_states if no state satisfies the condition.
    int get_min_state_above_energy(double energy) const;

    // Extract local k-point indices, used for MPI parallelization
    std::vector<int> get_iks_local() const;

    //! Get the density matrix of a particular spin and kpoint
    ComplexMatrix get_dmat_cplx(int ispin, int ispinor_bra, int ispinor_ket, int ikpt) const;

    // Density matrix and green's function calculation, serial version
    ComplexMatrix get_dmat_cplx_R(int ispin, int ispinor_bra, int ispinor_ket,
                                  const std::vector<Vector3_Order<double>>& kfrac_list,
                                  const Vector3_Order<int>& R) const;
    std::map<Vector3_Order<int>, ComplexMatrix> get_dmat_cplx_Rs(
        int ispin, int ispinor_bra, int ispinor_ket,
        const std::vector<Vector3_Order<double>>& kfrac_list,
        const std::vector<Vector3_Order<int>>& Rs) const;

    //! Optional band mask: true includes a band; an empty mask includes all bands.
    //! A nonempty mask must have get_n_bands() entries, shared by all k-points and spins.
    ComplexMatrix get_gf_cplx_imagtime(int ispin, int ispinor_bra, int ispinor_ket, int ikpt,
                                       double tau, const std::vector<bool> &band_mask = {}) const;
    std::map<double, std::map<Vector3_Order<int>, ComplexMatrix>> get_gf_cplx_imagtimes_Rs(
        int ispin, int ispinor_bra, int ispinor_ket,
        const std::vector<Vector3_Order<double>>& kfrac_list, std::vector<double> imagtimes,
        const std::vector<Vector3_Order<int>>& Rs, const std::vector<bool> &band_mask = {}) const;
    std::map<double, std::map<Vector3_Order<int>, matrix>> get_gf_real_imagtimes_Rs(
        int ispin, int ispinor_bra, int ispinor_ket,
        const std::vector<Vector3_Order<double>>& kfrac_list, std::vector<double> imagtimes,
        const std::vector<Vector3_Order<int>>& Rs, const std::vector<bool> &band_mask = {}) const;

    // void allredue_wfc_isk();
};

bool can_restore_symmetry_kstar_meanfield(
    const SymmetryContext& ctx,
    const std::vector<SpeciesBasisLayout>& wfc_layouts,
    const MeanField& mf,
    const std::vector<Vector3_Order<double>>& kfrac_list,
    const std::map<atom_t, size_t>& atom_nw);

symmetry_kstar_member_kfrac_targets_t build_symmetry_kstar_member_kfrac_targets(
    const SymmetryContext& ctx,
    const PeriodicBoundaryData& pbc);

ComplexMatrix get_symmetry_restored_dmat_cplx_R(
    const SymmetryContext& ctx,
    const std::vector<SpeciesBasisLayout>& wfc_layouts,
    const MeanField& mf,
    int ispin, int ispinor_bra, int ispinor_ket,
    const std::vector<Vector3_Order<double>>& kfrac_list,
    const Vector3_Order<int>& R,
    const std::map<atom_t, size_t>& atom_nw,
    const symmetry_kstar_member_kfrac_targets_t* member_kfrac_targets = nullptr,
    const symmetry_kstar_representative_indices_t* representative_k_indices = nullptr);

/*!
 * @brief Four-channel spinor variant of get_symmetry_restored_dmat_cplx_R.
 *
 * Returns the four spin blocks D^{ab}(R) (a = bra, b = ket,
 * channel-outermost convention C8) restored from the IBZ representatives.
 * Each star member is resolved through
 * `member.action_id -> ctx.kspace_actions -> ctx.spin_operations[canonical]`
 * into (spatial_id, U_s, eta), exactly as in the spinor Green's-function
 * restore. The orbital transform matrix A is built once per member (gauge
 * phases excluded) and reused for all four source blocks; the SU(2) mixing
 * and, for antiunitary members, the Theta remap
 * {conj(Y11), -conj(Y10), -conj(Y01), conj(Y00)} are applied by
 * transform_spinor_bilinear. The target-kpoint gauge phases are plain unitary
 * re-gauging factors applied after the kernel. The per-channel scalar restore
 * is NOT used here: it would miss the SU(2) block mixing and the channel
 * exchange of the Theta remap for antiunitary members.
 *
 * Requires mf.get_n_spinor() == 2. Missing (bra, ket) source blocks are
 * zero-filled (get_dmat_cplx returns a zero block for an absent wfc channel),
 * so a spin-mixing member cannot drop sparse blocks.
 */
SpinorBlocks4<ComplexMatrix> get_symmetry_restored_dmat_cplx_R_spinor(
    const SymmetryContext& ctx, const std::vector<SpeciesBasisLayout>& wfc_layouts,
    const MeanField& mf, int ispin, const std::vector<Vector3_Order<double>>& kfrac_list,
    const Vector3_Order<int>& R, const std::map<atom_t, size_t>& atom_nw,
    const symmetry_kstar_member_kfrac_targets_t* member_kfrac_targets = nullptr,
    const symmetry_kstar_representative_indices_t* representative_k_indices = nullptr);

std::map<double, std::map<Vector3_Order<int>, ComplexMatrix>>
get_symmetry_restored_gf_cplx_imagtimes_Rs(
    const SymmetryContext& ctx,
    const std::vector<SpeciesBasisLayout>& wfc_layouts,
    const MeanField& mf,
    int ispin, int ispinor_bra, int ispinor_ket,
    const std::vector<Vector3_Order<double>>& kfrac_list,
    const std::vector<double>& imagtimes,
    const std::vector<Vector3_Order<int>>& Rs,
    const std::map<atom_t, size_t>& atom_nw,
    const std::vector<bool> &band_mask = {},
    const symmetry_kstar_member_kfrac_targets_t* member_kfrac_targets = nullptr,
    const symmetry_kstar_representative_indices_t* representative_k_indices = nullptr);

/*!
 * @brief Validate that a Green's-function band cutoff does not slice through a
 * degenerate band multiplet on the k-point set used as symmetry-restore source.
 *
 * A symmetry restore mixes the selected bands across each k-star; when any
 * band-mask boundary separates two (nearly) degenerate eigenvalues the selected
 * space is not closed under the star operations and the restored GF would
 * depend on the arbitrary gauge inside the degenerate subspace. Throws
 * LIBRPA_RUNTIME_ERROR naming the k-point, the band indices and the gap when
 * the gap at any mask boundary is smaller than degen_tol.
 *
 * No-op when band_mask is empty, which denotes the full band space.
 */
void validate_kstar_band_cutoff_closure(
    const SymmetryContext& ctx,
    const MeanField& mf,
    const std::vector<bool>& band_mask,
    double degen_tol = 1e-8);

/*!
 * @brief Four-channel spinor variant of get_symmetry_restored_gf_cplx_imagtimes_Rs.
 *
 * Returns, for each tau and R, the four spin blocks G^{ab} (a = bra, b = ket,
 * channel-outermost convention C8) restored from the IBZ representatives.
 * Each star member is resolved through
 * `member.action_id -> ctx.kspace_actions -> ctx.spin_operations[canonical]`
 * into (spatial_id, U_s, eta). The orbital transform matrix A is built once
 * per member (gauge phases excluded) and reused for all four source blocks;
 * the SU(2) mixing and, for antiunitary members, the Theta remap
 * {conj(Y11), -conj(Y10), -conj(Y01), conj(Y00)} are applied by
 * transform_spinor_bilinear. The target-kpoint gauge phases are plain unitary
 * re-gauging factors applied after the kernel, so antiunitary members match
 * the scalar restore convention. Transforms act in the tau domain point by
 * point (tau is real, z* = tau; report section 6.8).
 *
 * Requires mf.get_n_spinor() == 2. Missing (bra, ket) source blocks are
 * zero-filled. Non-empty band_mask is guarded by
 * validate_kstar_band_cutoff_closure.
 */
std::map<double, std::map<Vector3_Order<int>, SpinorBlocks4<ComplexMatrix>>>
get_symmetry_restored_gf_cplx_imagtimes_Rs_spinor(
    const SymmetryContext& ctx,
    const std::vector<SpeciesBasisLayout>& wfc_layouts,
    const MeanField& mf,
    int ispin,
    const std::vector<Vector3_Order<double>>& kfrac_list,
    const std::vector<double>& imagtimes,
    const std::vector<Vector3_Order<int>>& Rs,
    const std::map<atom_t, size_t>& atom_nw,
    const std::vector<bool>& band_mask = {},
    const symmetry_kstar_member_kfrac_targets_t* member_kfrac_targets = nullptr,
    const symmetry_kstar_representative_indices_t* representative_k_indices = nullptr);

/*!
 * @brief Extract one (bra, ket) channel out of a four-channel spinor GF map,
 * moving the blocks out of the input.
 */
std::map<double, std::map<Vector3_Order<int>, ComplexMatrix>>
extract_spinor_gf_block(
    std::map<double, std::map<Vector3_Order<int>, SpinorBlocks4<ComplexMatrix>>>&& gf_spinor,
    int ispinor_bra, int ispinor_ket);

}

#endif // ! MEANFIELD_H
