#include "meanfield.h"
#include "../mpi/base_mpi.h"

namespace librpa_int
{

// Density matrix and green's function calculation, MPI parallel (over k-points) version

// User should ensure that wave function at a certain k-point exists only on one MPI process
std::map<Vector3_Order<int>, ComplexMatrix> get_dmat_cplx_Rs_kpara(
    int ispin, const MeanField &mf, const std::vector<Vector3_Order<double>>& kfrac_list,
    const std::vector<Vector3_Order<int>>& Rs, const MpiCommHandler &comm_h);

}
