#include <algorithm>
#include "constants.h"
#include "coulmat.h"
#include "pbc.h"

atpair_R_mat_t
FT_Vq(const atpair_k_cplx_mat_t &coulmat_k, vector<Vector3_Order<int>> Rlist, bool return_ordered_atom_pair)
{
    atpair_R_mat_t coulmat_R;

    for (auto R: Rlist)
    {
        auto iteR = std::find(Rlist.cbegin(), Rlist.cend(), R);
        auto iR = std::distance(Rlist.cbegin(), iteR);
        for (const auto &Mu_NuqV: coulmat_k)
        {
            const auto Mu = Mu_NuqV.first;
            const int n_mu = atom_mu[Mu];
            for (const auto &Nu_qV: Mu_NuqV.second)
            {
                const auto Nu = Nu_qV.first;
                const int n_nu = atom_mu[Nu];
                coulmat_R[Mu][Nu][R] = make_shared<matrix>();
                // a temporary complex matrix to save the transformed matrix
                ComplexMatrix VR_cplx(n_mu, n_nu);
                for (const auto &q_V: Nu_qV.second)
                {
                    auto q = q_V.first;
                    for (auto q_bz: map_irk_ks[q])
                    {
                        double ang = - q_bz * (R * latvec) * TWO_PI;
                        complex<double> kphase = complex<double>(cos(ang), sin(ang));
                        // FIXME: currently support inverse symmetry only
                        if (q_bz == q)
                        {
                            // cout << "Direct:  " << q_bz << " => " << q << ", phase = " << kphase << endl;
                            VR_cplx += (*q_V.second) * kphase;
                        }
                        else
                        {
                            // cout << "Inverse: " << q_bz << " => " << q << ", phase = " << kphase << endl;
                            VR_cplx += conj(*q_V.second) * kphase;
                        }
                    }
                    // minyez debug: check hermicity of Vq
                    // if (iR == 0)
                    // {
                    //     int iq = std::distance(klist.begin(), std::find(klist.begin(), klist.end(), q));
                    //     sprintf(fn, "Vq_Mu_%zu_Nu_%zu_iq_%d.mtx", Mu, Nu, iq);
                    //     print_complex_matrix_mm(*q_V.second, fn);
                    // }
                    // end minyez debug
                }
                *coulmat_R[Mu][Nu][R] = VR_cplx.real();
                // debug print
                // sprintf(fn, "VR_cplx_Mu_%zu_Nu_%zu_iR_%zu.mtx", Mu, Nu, iR);
                // print_complex_matrix_mm(VR_cplx, fn);
                // sprintf(fn, "VR_Mu_%zu_Nu_%zu_iR_%zu.mtx", Mu, Nu, iR);
                // print_matrix_mm(*coulmat_R[Mu][Nu][R], fn);

                // when ordered atom pair is requested, check whether it is available in the original map
                if (return_ordered_atom_pair && Mu != Nu && (coulmat_k.count(Nu) == 0 || coulmat_k.at(Nu).count(Mu) == 0))
                {
                    coulmat_R[Nu][Mu][R] = make_shared<matrix>();
                    ComplexMatrix VR_cplx(n_nu, n_mu);
                    for (const auto &q_V: Nu_qV.second)
                    {
                        auto q = q_V.first;
                        for (auto q_bz: map_irk_ks[q])
                        {
                            double ang = - q_bz * (R * latvec) * TWO_PI;
                            complex<double> kphase = complex<double>(cos(ang), sin(ang));
                            if (q_bz == q)
                            {
                                VR_cplx += transpose(*q_V.second, true) * kphase;
                            }
                            else
                            {
                                VR_cplx += transpose(*q_V.second, false) * kphase;
                            }
                        }
                    }
                    *coulmat_R[Nu][Mu][R] = VR_cplx.real();
                }
            }
        }
    }
    // myz debug: check the imaginary part of the coulomb matrix
    // char fn[80];
    /* for (const auto & Mu_NuRV: VR) */
    /* { */
    /*     auto Mu = Mu_NuRV.first; */
    /*     const int n_mu = atom_mu[Mu]; */
    /*     for (const auto & Nu_RV: Mu_NuRV.second) */
    /*     { */
    /*         auto Nu = Nu_RV.first; */
    /*         const int n_nu = atom_mu[Nu]; */
    /*         for (const auto & R_V: Nu_RV.second) */
    /*         { */
    /*             auto R = R_V.first; */
    /*             auto &V = R_V.second; */
    /*             auto iteR = std::find(Rlist.cbegin(), Rlist.cend(), R); */
    /*             auto iR = std::distance(Rlist.cbegin(), iteR); */
    /*             sprintf(fn, "VR_Mu_%zu_Nu_%zu_iR_%zu.mtx", Mu, Nu, iR); */
    /*             print_complex_matrix_mm(*V, fn); */
    /*         } */
    /*     } */
    /* } */
    return coulmat_R;
}

