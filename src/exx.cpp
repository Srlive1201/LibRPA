#include "params.h"
#include "constants.h"
#include "exx.h"
#include "input.h"
#include "lapack_connector.h"

void Exx::build_density_matrix_R(const MeanField &mf, const vector<Vector3_Order<double>> &klist, const Vector3_Order<int> &R)
{
    const auto nkpts = mf.get_n_kpoints();
    const auto nspins = mf.get_n_spins();
    const auto nbands = mf.get_n_bands();
    const auto naos = mf.get_n_aos();

    ComplexMatrix dmat_cplx(naos, naos);
    for (int is = 0; is != nspins; is++)
    {
        dmat_cplx.zero_out();
        for (int ik = 0; ik != nkpts; ik++)
        {
            double ang = - klist[ik] * (R * latvec) * TWO_PI;
            complex<double> kphase = complex<double>(cos(ang), sin(ang));
            dmat_cplx += kphase * mf.get_dmat_cplx(is, ik);
        }
        for ( int I = 0; I != natom; I++)
        {
            const auto I_num = atom_nw[I];
            for (int J = 0; J != natom; J++)
            {
                const size_t J_num = atom_nw[J];
                ComplexMatrix dmat_cplx_IJR(I_num, J_num);
                for (size_t i = 0; i != I_num; i++)
                {
                    size_t i_glo = atom_iw_loc2glo(I, i);
                    for (size_t j = 0; j != J_num; j++)
                    {
                        size_t j_glo = atom_iw_loc2glo(J, j);
                        dmat_cplx_IJR(i, j) = dmat_cplx(i_glo, j_glo);
                    }
                }
                if (dmat_cplx_IJR.get_max_imag())
                    printf("Warning: converting complex-valued Dmat to real, IJR %d %d (%d, %d, %d)", I, J, R.x, R.y, R.z);
                dmat[is][I][J][R] = std::make_shared<matrix>();
                *dmat[is][I][J][R] = dmat_cplx_IJR.real();
            }
        }
    }
}

void Exx::build_exx_orbital_energy(const MeanField &mf, const vector<Vector3_Order<double>> &klist, const atpair_R_mat_t &LRI_Cs, const vector<Vector3_Order<int>> &Rlist, const Vector3_Order<int> &R_period, const atpair_R_cplx_mat_t &coul_mat)
{

}
