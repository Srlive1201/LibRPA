#include "meanfield_mpi.h"

#include <algorithm>

#include "../math/matrix_m.h"
#include "../math/lapack_connector.h"
#include "../utils/constants.h"
#include "mpi_proto.h"

namespace librpa_int
{

std::map<Vector3_Order<int>, ComplexMatrix> get_dmat_cplx_Rs_kpara(
    int ispin, const MeanField &mf, const std::vector<Vector3_Order<double>>& kfrac_list,
    const std::vector<Vector3_Order<int>>& Rs, const MpiCommHandler &comm_h)
{
    std::map<Vector3_Order<int>, ComplexMatrix> dmat_local;
    const auto iks_local = mf.get_iks_local();
    const int nk_local = iks_local.size();
    const int n_aos = mf.get_n_aos();

    // Collect Rs requested by each process
    std::vector<int> n_Rs_all(comm_h.nprocs, 0);
    n_Rs_all[comm_h.myid] = Rs.size();
    MPI_Allreduce(MPI_IN_PLACE, n_Rs_all.data(), comm_h.nprocs, mpi_datatype<int>::value, MPI_SUM, comm_h.comm);
    const int nR_max = *std::max_element(n_Rs_all.cbegin(), n_Rs_all.cend());
    std::vector<int> Rs_all(comm_h.nprocs * nR_max * 3, 0);
    for (int iR = 0; iR < n_Rs_all[comm_h.myid]; iR++)
    {
        Rs_all[comm_h.myid * nR_max * 3 + iR * 3] = Rs[iR].x;
        Rs_all[comm_h.myid * nR_max * 3 + iR * 3 + 1] = Rs[iR].y;
        Rs_all[comm_h.myid * nR_max * 3 + iR * 3 + 2] = Rs[iR].z;
    }
    MPI_Allreduce(MPI_IN_PLACE, Rs_all.data(), comm_h.nprocs * nR_max * 3, mpi_datatype<int>::value, MPI_SUM, comm_h.comm);

    const size_t size = n_aos * n_aos;

    Matz kmat(size, nk_local, MAJOR::COL);
    Matz transmat(nk_local, nR_max, MAJOR::COL);
    Matz rmat(size, nR_max, MAJOR::COL);

    // NOTE: Support dimension (number of basis) up to 46300 (Gamma-only) or 1930 (k 8x8x8) due to range of int of axpy.
    //       Need division into batches for larger system, particularly for memory consideration
    for (int ik_local = 0; ik_local < nk_local; ik_local++)
    {
        const auto dmat_k = mf.get_dmat_cplx(ispin, iks_local[ik_local]);
        memcpy(kmat.ptr() + size * ik_local, dmat_k.c, size * sizeof(Matz::type));
    }

    for (int pid = 0; pid < comm_h.nprocs; pid++)
    {
        const auto nR_this = n_Rs_all[pid];
        rmat = 0.0;
        if (nR_this < 1) continue;
        for (int iR = 0; iR <nR_this; iR++)
        {
            const auto R_this = Rs_all.data() + pid * nR_max * 3 + iR * 3;
            for (int ik = 0; ik < nk_local; ik++)
            {
                const auto &kf = kfrac_list[iks_local[ik]];
                auto ang = - (kf.x * R_this[0] + kf.y * R_this[1] + kf.z * R_this[0]) * TWO_PI;
                transmat(ik, iR) = cplxdb{cos(ang), sin(ang)};
            }
        }
        if (nk_local > 0)
            LapackConnector::gemm_f('N', 'N', size, nR_this, nk_local, 1.0,
                                    kmat.ptr(), size, transmat.ptr(), nk_local, 0.0, rmat.ptr(), size);
        const size_t count = nk_local > 0? size * nR_this : 0;
        MPI_Reduce(MPI_IN_PLACE, rmat.ptr(), count, mpi_datatype<cplxdb>::value, MPI_SUM, pid, comm_h.comm);
        if (comm_h.myid == pid)
        {
            for (int iR = 0; iR <nR_this; iR++)
            {
                const auto R_this = Rs_all.data() + pid * nR_max * 3 + iR * 3;
                ComplexMatrix m(n_aos, n_aos);
                memcpy(m.c, rmat.ptr() + size * iR, size * sizeof(cplxdb));
                Vector3_Order<int> R{R_this[0], R_this[1], R_this[2]};
                dmat_local.emplace(R, m);
            }
        }
        comm_h.barrier(); // May try out non-blocking later
    }

    return dmat_local;
}

}
