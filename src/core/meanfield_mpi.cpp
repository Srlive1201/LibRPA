#include "meanfield_mpi.h"

#include <algorithm>
#include <complex>

#include "../io/stl_io_helper.h"
#include "../math/lapack_connector.h"
#include "../math/matrix_m.h"
#include "../utils/constants.h"

namespace librpa_int
{

static void collect_Rs(const std::vector<Vector3_Order<int>> &Rs, std::vector<int> &n_Rs_all,
                       std::vector<int> &Rs_all, int &nR_max, const MpiCommHandler &comm_h)
{
    n_Rs_all.resize(comm_h.nprocs);
    for (int pid = 0; pid < comm_h.nprocs; pid++) n_Rs_all[pid] = 0;
    const int n_Rs_this = Rs.size();
    n_Rs_all[comm_h.myid] = n_Rs_this;
    // global::ofs_myid << global::myid_global << " " << global::size_global << std::endl;
    // comm_h.barrier();
    MPI_Allreduce(MPI_IN_PLACE, n_Rs_all.data(), comm_h.nprocs, mpi_datatype<int>::value, MPI_SUM, comm_h.comm);
    nR_max = *std::max_element(n_Rs_all.cbegin(), n_Rs_all.cend());
    Rs_all.resize(3 * nR_max * comm_h.nprocs);
    for (int iR = 0; iR < n_Rs_this; iR++)
    {
        Rs_all[comm_h.myid * nR_max * 3 + iR * 3] = Rs[iR].x;
        Rs_all[comm_h.myid * nR_max * 3 + iR * 3 + 1] = Rs[iR].y;
        Rs_all[comm_h.myid * nR_max * 3 + iR * 3 + 2] = Rs[iR].z;
    }
    MPI_Allreduce(MPI_IN_PLACE, Rs_all.data(), comm_h.nprocs * nR_max * 3, mpi_datatype<int>::value, MPI_SUM, comm_h.comm);
}

std::map<Vector3_Order<int>, ComplexMatrix> get_dmat_cplx_Rs_kpara(
    int ispin, const MeanField &mf, const std::vector<Vector3_Order<double>>& kfrac_list,
    const std::vector<Vector3_Order<int>>& Rs, const MpiCommHandler &comm_h)
{
    std::map<Vector3_Order<int>, ComplexMatrix> dmat_local;

    // Collect Rs requested by each process
    std::vector<int> n_Rs_all, Rs_all;
    int nR_max;
    collect_Rs(Rs, n_Rs_all, Rs_all, nR_max, comm_h);
    // global::ofs_myid << "get_dmat_cplx_Rs_kpara nRs_all " << n_Rs_all << std::endl;

    const int n_aos = mf.get_n_aos();
    const size_t size = n_aos * n_aos;

    const auto iks_local = mf.get_iks_local();
    // global::ofs_myid << "iks_local " << iks_local << std::endl;
    const int nk_local = iks_local.size();
    // Check if there is duplicate k-point data
    int nk_local_sum;
    MPI_Allreduce(&nk_local, &nk_local_sum, 1, MPI_INT, MPI_SUM, comm_h.comm);
    if (nk_local_sum > mf.get_n_kpoints())
        throw LIBRPA_RUNTIME_ERROR("found duplicated k-point eigenvectors data");
    else if (nk_local_sum < mf.get_n_kpoints())
        throw LIBRPA_RUNTIME_ERROR("missing k-point eigenvectors data");

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
        // NOTE: MPI_Reduce requires all processes use the same sendcount, but rmat is simply zero where nk_local == 0.
        const size_t count = size * nR_this;
        // global::ofs_myid << "get_dmat_cplx_Rs_kpara pid " << pid << " nR_this " << nR_this << " count " << count << " matsize " << rmat.size() << std::endl;

        if (nR_this < 1) continue;
        #pragma omp parallel for collapse(2) schedule(dynamic)
        for (int iR = 0; iR < nR_this; iR++)
        {
            for (int ik = 0; ik < nk_local; ik++)
            {
                const int index = pid * nR_max * 3 + iR * 3;
                const auto &kf = kfrac_list[iks_local[ik]];
                auto ang = - (kf.x * Rs_all[index] + kf.y * Rs_all[index+1] + kf.z * Rs_all[index+2]) * TWO_PI;
                transmat(ik, iR) = cplxdb{cos(ang), sin(ang)};
            }
        }
        // global::ofs_myid << "transmat pid" << std::endl;
        // global::ofs_myid << transmat << std::endl;
        rmat = 0.0;
        if (nk_local > 0)
        {
            LapackConnector::gemm_f('N', 'N', size, nR_this, nk_local, 1.0,
                                    kmat.ptr(), size, transmat.ptr(), nk_local, 0.0, rmat.ptr(), size);
        }
        // global::ofs_myid << rmat << std::endl;
        if (comm_h.myid == pid)
        {
            MPI_Reduce(MPI_IN_PLACE, rmat.ptr(), count, mpi_datatype<Matz::type>::value, MPI_SUM, pid, comm_h.comm);
            for (int iR = 0; iR < nR_this; iR++)
            {
                const int index = pid * nR_max * 3 + iR * 3;
                ComplexMatrix m(n_aos, n_aos);
                memcpy(m.c, rmat.ptr() + size * iR, size * sizeof(Matz::type));
                Vector3_Order<int> R{Rs_all[index], Rs_all[index+1], Rs_all[index+2]};
                // global::ofs_myid << "iR " << iR << " R " << R << std::endl;
                dmat_local.emplace(R, std::move(m));
            }
        }
        else
        {
            MPI_Reduce(rmat.ptr(), rmat.ptr(), count, mpi_datatype<cplxdb>::value, MPI_SUM, pid, comm_h.comm);
        }
        comm_h.barrier(); // May try out non-blocking later
    }

    return dmat_local;
}

std::map<double, std::map<Vector3_Order<int>, ComplexMatrix>> get_gf_cplx_imagtimes_Rs_kpara(
    int ispin, const MeanField &mf, const std::vector<Vector3_Order<double>> &kfrac_list, std::vector<double> imagtimes,
    const std::vector<Vector3_Order<int>> &Rs, const MpiCommHandler &comm_h)
{
    std::map<double, std::map<Vector3_Order<int>, ComplexMatrix>> gf;

    // Collect Rs requested by each process
    std::vector<int> n_Rs_all, Rs_all;
    int nR_max;
    collect_Rs(Rs, n_Rs_all, Rs_all, nR_max, comm_h);
    // global::ofs_myid << "get_gf_cplx_imagtimes_Rs_kpara nRs_all " << n_Rs_all << std::endl;

    const int n_aos = mf.get_n_aos();
    const size_t size = n_aos * n_aos;

    const auto iks_local = mf.get_iks_local();
    // global::ofs_myid << "iks_local " << iks_local << std::endl;
    // Check if there is duplicate k-point data
    const int nk_local = iks_local.size();
    int nk_local_sum;
    MPI_Allreduce(&nk_local, &nk_local_sum, 1, MPI_INT, MPI_SUM, comm_h.comm);
    if (nk_local_sum > mf.get_n_kpoints())
        throw LIBRPA_RUNTIME_ERROR("found duplicated k-point eigenvectors data");
    else if (nk_local_sum < mf.get_n_kpoints())
        throw LIBRPA_RUNTIME_ERROR("missing k-point eigenvectors data");

    Matz kmat(size, nk_local, MAJOR::COL);
    Matz transmat(nk_local, nR_max, MAJOR::COL);
    Matz rmat(size, nR_max, MAJOR::COL);

    for (auto tau: imagtimes)
    {
        std::map<Vector3_Order<int>, ComplexMatrix> gf_tau;
        // TODO: this part is the same as denstiy matrix calculation, so it may be extracted to a common function
        for (int ik_local = 0; ik_local < nk_local; ik_local++)
        {
            const auto dmat_k = mf.get_gf_cplx_imagtime(ispin, iks_local[ik_local], tau);
            memcpy(kmat.ptr() + size * ik_local, dmat_k.c, size * sizeof(Matz::type));
        }

        for (int pid = 0; pid < comm_h.nprocs; pid++)
        {
            const auto nR_this = n_Rs_all[pid];
            // NOTE: MPI_Reduce requires all processes use the same sendcount, but rmat is simply zero where nk_local == 0.
            const size_t count = size * nR_this;
            // global::ofs_myid << "get_dmat_cplx_Rs_kpara pid " << pid << " nR_this " << nR_this << " count " << count << " matsize " << rmat.size() << std::endl;

            if (nR_this < 1) continue;
            #pragma omp parallel for collapse(2) schedule(dynamic)
            for (int iR = 0; iR <nR_this; iR++)
            {
                for (int ik = 0; ik < nk_local; ik++)
                {
                    const auto R_this = Rs_all.data() + pid * nR_max * 3 + iR * 3;
                    const auto &kf = kfrac_list[iks_local[ik]];
                    auto ang = - (kf.x * R_this[0] + kf.y * R_this[1] + kf.z * R_this[0]) * TWO_PI;
                    transmat(ik, iR) = cplxdb{cos(ang), sin(ang)};
                }
            }
            rmat = 0.0;
            if (nk_local > 0)
            {
                LapackConnector::gemm_f('N', 'N', size, nR_this, nk_local, 1.0,
                                        kmat.ptr(), size, transmat.ptr(), nk_local, 0.0, rmat.ptr(), size);
            }
            // global::ofs_myid << rmat << std::endl;
            if (comm_h.myid == pid)
            {
                MPI_Reduce(MPI_IN_PLACE, rmat.ptr(), count, mpi_datatype<cplxdb>::value, MPI_SUM, pid, comm_h.comm);
                for (int iR = 0; iR < nR_this; iR++)
                {
                    const auto R_this = Rs_all.data() + pid * nR_max * 3 + iR * 3;
                    ComplexMatrix m(n_aos, n_aos);
                    memcpy(m.c, rmat.ptr() + size * iR, size * sizeof(cplxdb));
                    Vector3_Order<int> R{R_this[0], R_this[1], R_this[2]};
                    gf_tau.emplace(R, std::move(m));
                }
            }
            else
            {
                MPI_Reduce(rmat.ptr(), rmat.ptr(), count, mpi_datatype<cplxdb>::value, MPI_SUM, pid, comm_h.comm);
            }
            comm_h.barrier(); // May try out non-blocking later
        }
        gf.emplace(tau, std::move(gf_tau));
    }

    return gf;
}

}
