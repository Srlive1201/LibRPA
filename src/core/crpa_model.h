#pragma once

#include <array>
#include <utility>

#include "crpa_workflow.h"

namespace librpa_int
{
//! Owning electronic-state snapshot. All complex matrices are row major.
struct CrpaModelK
{
    std::vector<double> energies, occupations;
    ComplexMatrix coefficients;  //!< C[component-major AO][band], not stored wfc=C^T.
    ComplexMatrix parent_phi, output_phi;
    ComplexMatrix parent_projector, output_projector;  //!< P=Phi^dagger S C.
    ComplexMatrix parent_h, output_h;                  //!< P diag(E) P^dagger, in Hartree.
    ComplexMatrix parent_density, output_density;      //!< P diag(f) P^dagger.
    double parent_capture_min = 1.0;
};

struct CrpaModelSnapshot
{
    int n_spins = 0, n_spinor = 0, n_aos = 0, n_bands = 0;
    double fermi = 0;
    std::vector<ComplexMatrix> overlap;
    std::vector<std::vector<CrpaModelK>> spin_k;
    std::vector<SiteOrbitalGroup> parent_sites, output_sites;
    std::vector<int> output_columns;
};

//! No eigenproblem or response calculation: retain the actual contraction frame.
CrpaModelSnapshot snapshot_crpa_model(CrpaPreparedModel& model, const MeanField& meanfield);
ComplexMatrix project_crpa_one_body(const ComplexMatrix& projector,
                                    const std::vector<double>& diagonal);
ComplexMatrix project_crpa_orbitals(const ComplexMatrix& coefficients,
                                    const ComplexMatrix& spatial_overlap, const ComplexMatrix& phi,
                                    int n_spinor);
ComplexMatrix rotate_crpa_projector(const ComplexMatrix& projector, const ComplexMatrix& rotation);
//! Phi_new=Phi_old R; both R matrices must have orthonormal columns.
ComplexMatrix rotate_crpa_tensor(const ComplexMatrix& tensor, const ComplexMatrix& left_rotation,
                                 const ComplexMatrix& right_rotation);
ComplexMatrix crpa_tensor_to_fermion(const ComplexMatrix& tensor, int n_orbitals);
//! Ordered-column means U_mm,mm, U_mm,nn and U_mn,nm; off-diagonal means vanish for M=1.
std::array<std::complex<double>, 3> crpa_tensor_averages(const ComplexMatrix& tensor,
                                                         int n_orbitals);
std::pair<ComplexMatrix, ComplexMatrix> average_crpa_local_model(
    const CrpaModelSnapshot& model, bool parent, int spin, int site,
    const std::vector<double>& k_weights);
ComplexMatrix fourier_crpa_model_hopping(const CrpaModelSnapshot& model, bool parent, int spin,
                                         const std::vector<Vector3_Order<double>>& k_fractional,
                                         const std::vector<double>& k_weights,
                                         const Vector3_Order<int>& lattice_R);
}  // namespace librpa_int
