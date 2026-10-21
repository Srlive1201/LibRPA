#pragma once

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
//! Phi_new=Phi_old R; both R matrices must have orthonormal columns.
ComplexMatrix rotate_crpa_tensor(const ComplexMatrix& tensor, const ComplexMatrix& left_rotation,
                                 const ComplexMatrix& right_rotation);
}  // namespace librpa_int
