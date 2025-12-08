#include "../core/meanfield_mpi.h"

#include <iostream>
#include <iterator>

#include "../mpi/global_mpi.h"
#include "../io/stl_io_helper.h"
#include "testutils.h"

using namespace std;
using namespace librpa_int;
using namespace librpa_int::global;

static MeanField init_mf_pbc(int nk, int nb, int nocc, double gap,
                             std::vector<Vector3_Order<double>> &kfrac_list,
                             std::vector<Vector3_Order<int>> &Rs)
{
    assert(nk > 0 && nb > 0);
    const int nspin = 1;
    MeanField mf(nspin, nk, nb, nb);
    const double efermi = 0.0;
    mf.get_efermi() = efermi;
    gap = std::abs(gap);
    for (int ik = 0; ik < nk; ik++)
    {
        for (int ib = 0; ib < nocc; ib++)
        {
            mf.get_weight().at(0)(ik, ib) = 2.0 / nspin / nk;
            mf.get_eigenvals().at(0)(ik, ib) = efermi - 0.5 * gap;
        }
        for (int ib = nocc; ib < nb; ib++)
        {
            mf.get_eigenvals().at(0)(ik, ib) = efermi + 0.5 * gap;
        }
    }

    kfrac_list.clear();
    Rs.clear();
    for (int ik = 0; ik < nk; ik++)
    {
        Vector3_Order<double> kf{double(ik) / nk, 0, 0};
        kfrac_list.emplace_back(kf);
        if (ik % size_global == myid_global)
        {
            ComplexMatrix mat(nb, nb);
            // Every eigenvector matrix is identity
            // Only (0,0,0) has non-zero density matrix. Other R-point are cancelled out.
            mat.set_as_identity_matrix();
            mf.get_eigenvectors()[0][ik] = std::move(mat);
        }
    }

    for (int iR = 0; iR < nk; iR++)
    {
        Vector3_Order<int> R{iR, 0, 0};
        if (iR % size_global == myid_global) Rs.emplace_back(R);
    }

    return mf;
}

static void test_dmat_cplx_Rs_kpara(int nk, int nb, int nocc)
{
    // ofs_myid << kfrac_list << endl;
    std::vector<Vector3_Order<double>> kfrac_list;
    std::vector<Vector3_Order<int>> Rs;
    const auto mf = init_mf_pbc(nk, nb, nocc, 1.0, kfrac_list, Rs);

    const auto dm_Rs = get_dmat_cplx_Rs_kpara(0, mf, kfrac_list, Rs, mpi_comm_global_h);
    assert(dm_Rs.size() == Rs.size());
    for (const auto &[R, rmat]: dm_Rs)
    {
        if (R.x == 0 && R.y == 0 && R.z == 0)
        {
            for (int ib = 0; ib < nocc; ib++)
            {
                assert(fequal(rmat(ib, ib), {1.0, 0.0}));
            }
            for (int ib = nocc; ib < nb; ib++)
            {
                assert(fequal(rmat(ib, ib), {0.0, 0.0}));
            }
        }
        else
        {
            for (int ib = 0; ib < nb; ib++)
            {
                assert(fequal(rmat(ib, ib), {0.0, 0.0}));
            }
        }
    }
}

static void test_gf_cplx_Rs_kpara(int nk, int nb, int nocc, double gap, double tau)
{
    std::vector<Vector3_Order<double>> kfrac_list;
    std::vector<Vector3_Order<int>> Rs;
    const auto mf = init_mf_pbc(nk, nb, nocc, gap, kfrac_list, Rs);

    const auto gf_Rs = get_gf_cplx_imagtimes_Rs_kpara(0, mf, kfrac_list, {tau}, Rs, mpi_comm_global_h);
    assert(gf_Rs.size() == 1);
    assert(gf_Rs.at(tau).size() == Rs.size());
    const double coeff = std::exp(- 0.5 * std::abs(gap * tau)) * (tau > 0? 1.0: -1.0);
    for (const auto &[tau, trmat]: gf_Rs)
    {
        for (const auto &[R, rmat]: trmat)
        {
            // ofs_myid << "tau " << tau << " R " << R << endl;
            // print_complex_matrix("", rmat, ofs_myid, true);
            if (R.x == 0 && R.y == 0 && R.z == 0)
            {
                for (int ib = 0; ib < nocc; ib++)
                {
                    if (tau > 0) assert(fequal(rmat(ib, ib), {0.0, 0.0}));
                    else assert(fequal(rmat(ib, ib), {coeff, 0.0}));
                }
                for (int ib = nocc; ib < nb; ib++)
                {
                    if (tau > 0) assert(fequal(rmat(ib, ib), {coeff, 0.0}));
                    else assert(fequal(rmat(ib, ib), {0.0, 0.0}));
                }
            }
            else
            {
                for (int ib = 0; ib < nb; ib++)
                {
                    assert(fequal(rmat(ib, ib), {0.0, 0.0}));
                }
            }
        }
    }
}
int main (int argc, char *argv[])
{
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    init_global_mpi();
    init_global_io();

    if ( size_global != 4 )
        throw std::runtime_error("test imposes 4 MPI processes");

    test_dmat_cplx_Rs_kpara(3, 8, 2);
    test_dmat_cplx_Rs_kpara(15, 4, 2);

    test_gf_cplx_Rs_kpara(3, 7, 2, 1.0, 1.0);
    test_gf_cplx_Rs_kpara(15, 3, 2, 0.5, -1.0);

    finalize_global_io();
    finalize_global_mpi();
    MPI_Finalize();

    return 0;
}
