#pragma once

#include <string>

#include "meanfield.h"

#include "../math/matrix_m.h"
#include "../mpi/base_mpi.h"
#include "../mpi/kpoint_blacs_parallel_context.h"

namespace librpa_int
{

struct SymmetryExpandedMeanField
{
    MeanField meanfield;
    std::vector<int> source_k;
    std::vector<bool> antiunitary;
};

// Expand paired source-gauge IBZ states without a new eigensolve. Returned
// band-space operators copy their source block, conjugating antiunitary routes.
SymmetryExpandedMeanField expand_symmetry_meanfield(const MeanField &source,
                                                    const AtomicBasis &basis,
                                                    const PeriodicBoundaryData &pbc,
                                                    const SymmetryContext &context);

// Verify E, occupations and diagonal endpoint weights through magnetic sewing.
// Degenerate permutations/mixing are allowed only when they preserve the weights.
// Includes all active operations, including little-group stabilizers.
// Retained rectangular frames must span a symmetry-invariant subspace. With no
// endpoint weights, an IBZ state is temporarily expanded only for validation.
void validate_symmetry_band_weights(
    const MeanField &mf, const AtomicBasis &basis, const PeriodicBoundaryData &pbc,
    const SymmetryContext &symmetry,
    const std::vector<std::vector<std::vector<double>>> &weights = {});

//! Shared block-size cap for permanent k-BLACS WFC layouts and their GEMM workspaces.
constexpr int wfc_gemm_block_size_opt = 128;

//! Redistribute every locally owned k-point wave function from one k-BLACS
//! layout to another and replace MeanField::wfc with the destination blocks.
void redistribute_meanfield_eigvecs_kblacs(
    MeanField &mf, const KPointBlacsParallelContext &kblacs_ctxt,
    const ArrayDesc &desc_wfc_src, const ArrayDesc &desc_wfc_dst,
    const std::string &label);

// ==================================================================
// Density matrix and green's function calculation, MPI parallel (over k-points) version
// User should ensure that wave function at a certain k-point exists only on one MPI process

std::map<Vector3_Order<int>, ComplexMatrix> get_dmat_cplx_Rs_kpara(
    int ispin, int ispinor_bra, int ispinor_ket, const MeanField &mf,
    const std::vector<Vector3_Order<double>> &kfrac_list, const std::vector<Vector3_Order<int>> &Rs,
    const MpiCommHandler &comm_h);

std::map<Vector3_Order<int>, ComplexMatrix> get_dmat_cplx_Rs_kpara(
    int ispin, const MeanField &mf,
    const std::vector<Vector3_Order<double>> &kfrac_list, const std::vector<Vector3_Order<int>> &Rs,
    const MpiCommHandler &comm_h);

std::map<double, std::map<Vector3_Order<int>, ComplexMatrix>> get_gf_cplx_imagtimes_Rs_kpara(
    int ispin, int ispinor_bra, int ispinor_ket, const MeanField &mf,
    const std::vector<Vector3_Order<double>> &kfrac_list, std::vector<double> imagtimes,
    const std::vector<Vector3_Order<int>> &Rs, const MpiCommHandler &comm_h, const std::vector<bool> &band_mask = {});

std::map<double, std::map<Vector3_Order<int>, ComplexMatrix>> get_gf_cplx_imagtimes_Rs_kpara(
    int ispin, const MeanField &mf,
    const std::vector<Vector3_Order<double>> &kfrac_list, std::vector<double> imagtimes,
    const std::vector<Vector3_Order<int>> &Rs, const MpiCommHandler &comm_h, const std::vector<bool> &band_mask = {});

// ==================================================================
// Density matrix and green's function calculation, two-level parallel (over k-points and BLACS) version.
// User should ensure that wave function are distributed among processes in kblacs_ctxt,
// and there is no duplicated k-point data.

// The output complex matrices are column-major Matz, which is different from conventions of functions above.
// This is mainly to account for usage of column major in BLACS functions, and can help to perform 2D->IJ transformation later
// using IndexScheduler.

// Rs: real-space vectors requested on this process. Processes in the same comm_kpoint may request
//     different R lists while sharing the same local BLACS matrix region.
// desc_wfc: descriptor of eigenvectors source of the MeanField object
// desc_gf: descriptor of the resulting distributed density matrix
std::map<Vector3_Order<int>, Matz> get_dmat_cplx_Rs_kblacs_para(
    int ispin, int ispinor_bra, int ispinor_ket, const MeanField &mf,
    const std::vector<Vector3_Order<double>> &kfrac_list, const std::vector<Vector3_Order<int>> &Rs,
    const KPointBlacsParallelContext &kblacs_ctxt, const ArrayDesc &desc_wfc, const ArrayDesc &desc_dm);

std::map<Vector3_Order<int>, Matz> get_dmat_cplx_Rs_kblacs_para(
    int ispin, const MeanField &mf,
    const std::vector<Vector3_Order<double>> &kfrac_list, const std::vector<Vector3_Order<int>> &Rs,
    const KPointBlacsParallelContext &kblacs_ctxt, const ArrayDesc &desc_wfc, const ArrayDesc &desc_dm);

std::map<Vector3_Order<int>, Matz> get_symmetry_restored_dmat_cplx_Rs_kblacs_para(
    int ispin, int ispinor_bra, int ispinor_ket, const MeanField &mf,
    const std::vector<Vector3_Order<double>> &kfrac_list, const std::vector<Vector3_Order<int>> &Rs,
    const KPointBlacsParallelContext &kblacs_ctxt, const ArrayDesc &desc_wfc, const ArrayDesc &desc_dm,
    const SymmetryContext &symmetry_context, const PeriodicBoundaryData &pbc,
    const AtomicBasis &atbasis_wfc);

/*!
 * @brief Four-channel spinor variant of
 * get_symmetry_restored_dmat_cplx_Rs_kblacs_para.
 *
 * Requires mf.get_n_spinor() == 2 and a restorable k-star context (the caller
 * checks can_restore_symmetry_kstar_meanfield; this function throws otherwise).
 * Each local k-point builds the four (bra, ket) source blocks with the same
 * pgemm pattern as the scalar path (missing wfc channels are zero-filled with
 * a stand-in block so the collective stays intact); each star member gathers
 * and spatially rotates the four block sets (time reversal suppressed, gauge
 * phases excluded), then applies the SU(2) mixing and the Theta remap through
 * transform_spinor_bilinear per atom pair. The target-kpoint gauge phases are
 * applied after the kernel as plain unitary re-gauging factors. The four
 * channels are packed into one communication buffer, so the reduce rounds are
 * unchanged from the scalar path.
 */
std::map<Vector3_Order<int>, SpinorBlocks4<Matz>>
get_symmetry_restored_dmat_cplx_Rs_kblacs_para_spinor(
    int ispin, const MeanField &mf, const std::vector<Vector3_Order<double>> &kfrac_list,
    const std::vector<Vector3_Order<int>> &Rs, const KPointBlacsParallelContext &kblacs_ctxt,
    const ArrayDesc &desc_wfc, const ArrayDesc &desc_dm, const SymmetryContext &symmetry_context,
    const PeriodicBoundaryData &pbc, const AtomicBasis &atbasis_wfc);

// Similar to above, but for Green's function. Note that imagtimes must be the same across kblacs_ctxt.
std::map<double, std::map<Vector3_Order<int>, Matz>> get_gf_cplx_imagtimes_Rs_kblacs_para(
    int ispin, int ispinor_bra, int ispinor_ket, const MeanField &mf,
    const std::vector<Vector3_Order<double>> &kfrac_list, std::vector<double> imagtimes,
    const std::vector<Vector3_Order<int>> &Rs,
    const KPointBlacsParallelContext &kblacs_ctxt, const ArrayDesc &desc_wfc, const ArrayDesc &desc_dm,
    const std::vector<bool> &band_mask = {});

std::map<double, std::map<Vector3_Order<int>, Matz>>
get_symmetry_restored_gf_cplx_imagtimes_Rs_kblacs_para(
    int ispin, int ispinor_bra, int ispinor_ket, const MeanField &mf,
    const std::vector<Vector3_Order<double>> &kfrac_list, std::vector<double> imagtimes,
    const std::vector<Vector3_Order<int>> &Rs,
    const KPointBlacsParallelContext &kblacs_ctxt, const ArrayDesc &desc_wfc, const ArrayDesc &desc_dm,
    const SymmetryContext &symmetry_context, const PeriodicBoundaryData &pbc,
    const AtomicBasis &atbasis_wfc, const std::vector<bool> &band_mask = {});

/*!
 * @brief Four-channel spinor variant of
 * get_symmetry_restored_gf_cplx_imagtimes_Rs_kblacs_para.
 *
 * Requires mf.get_n_spinor() == 2 and a restorable k-star context (the caller
 * checks can_restore_symmetry_kstar_meanfield; this function throws otherwise).
 * Each local k-point builds the four (bra, ket) source blocks with the same
 * pgemm pattern as the scalar path; each star member gathers and spatially
 * rotates the four block sets (time reversal suppressed, gauge phases
 * excluded), then applies the SU(2) mixing and the Theta remap through
 * transform_spinor_bilinear per atom pair. The target-kpoint gauge phases
 * are applied after the kernel as plain unitary re-gauging factors. The four
 * channels are packed into one communication buffer, so the reduce rounds
 * are unchanged from the scalar path. Missing (bra, ket) source blocks are
 * zero-filled.
 */
std::map<double, std::map<Vector3_Order<int>, SpinorBlocks4<Matz>>>
get_symmetry_restored_gf_cplx_imagtimes_Rs_kblacs_para_spinor(
    int ispin, const MeanField &mf,
    const std::vector<Vector3_Order<double>> &kfrac_list, std::vector<double> imagtimes,
    const std::vector<Vector3_Order<int>> &Rs,
    const KPointBlacsParallelContext &kblacs_ctxt, const ArrayDesc &desc_wfc, const ArrayDesc &desc_dm,
    const SymmetryContext &symmetry_context, const PeriodicBoundaryData &pbc,
    const AtomicBasis &atbasis_wfc);

/*!
 * @brief Extract one (bra, ket) channel out of a four-channel spinor GF map
 * of BLACS-distributed matrices, moving the blocks out of the input.
 */
std::map<double, std::map<Vector3_Order<int>, Matz>>
extract_spinor_gf_block_kblacs(
    std::map<double, std::map<Vector3_Order<int>, SpinorBlocks4<Matz>>>&& gf_spinor,
    int ispinor_bra, int ispinor_ket);

std::map<double, std::map<Vector3_Order<int>, Matz>> get_gf_cplx_imagtimes_Rs_kblacs_para(
    int ispin, const MeanField &mf,
    const std::vector<Vector3_Order<double>> &kfrac_list, std::vector<double> imagtimes,
    const std::vector<Vector3_Order<int>> &Rs,
    const KPointBlacsParallelContext &kblacs_ctxt, const ArrayDesc &desc_wfc, const ArrayDesc &desc_dm,
    const std::vector<bool> &band_mask = {});

}  // namespace librpa_int
