/*! @file: chi0.cpp
 */
#include <iostream>
#include "parallel_mpi.h"
#include "chi0.h"
#include "pbc.h"
#include <omp.h>
#include "matrix.h"
#include "complexmatrix.h"
#include "lapack_connector.h"
#include "input.h"
#include "constants.h"

void Chi0::compute(const atpair_R_mat_t &LRI_Cs,
                   const vector<Vector3_Order<int>> &Rlist,
                   const Vector3_Order<int> &R_period,
                   const vector<atpair_t> &atpairs_ABF,
                   const map<int, Vector3_Order<double>> &qlist,
                   TFGrids::GRID_TYPES gt, bool use_space_time)
{
    gf_save = gf_discard = 0;
    // reset chi0_q in case the method was called before
    chi0_q.clear();

    // prepare the time-freq grids
    switch (gt)
    {
        case ( TFGrids::GRID_TYPES::Minimax):
        {
            double emin, emax;
            mf.get_E_min_max(emin, emax);
            tfg.generate_minimax(emin, emax);
            break;
        }
        case ( TFGrids::GRID_TYPES::EvenSpaced_TF ):
        {
            // WARN: only used to debug, adjust emin and interval manually
            tfg.generate_evenspaced_tf(0.0304601, 0.0, 0.0116073, 0.0);
            break;
        }
        default:
            throw invalid_argument("chi0 on requested grid is not implemented");
    }

    if (use_space_time && !tfg.has_time_grids())
        throw invalid_argument("current grid type does not have time grids");

    if (use_space_time)
    {
        // use space-time method
        for (int iR = 0; iR < Rlist.size(); iR++)
            for (int itau = 0; itau < tfg.size(); itau++)
            {
                build_gf_Rt(iR, Rlist[iR], itau, 'O');
                build_gf_Rt(iR, Rlist[iR], itau, 'V');
            }
        cout << " FINISH GREEN FUN" << endl;
        cout << " Green threshold: " << gf_R_threshold << endl;
        cout << " Green save num: " << gf_save << endl;
        cout << " Green discard num: " << gf_discard << endl;
        build_chi0_q_space_time(LRI_Cs, Rlist, R_period, atpairs_ABF, qlist);
    }
    else
    {
        // conventional method does not need to build Green's function explicitly
        build_chi0_q_conventional(LRI_Cs, Rlist, R_period, atpairs_ABF, qlist);
    }
}

void Chi0::build_gf_Rt(size_t iR, Vector3_Order<int> R, size_t itau, char ov)
{
    double tau;
    auto nkpts = mf.get_n_kpoints();
    auto nspins = mf.get_n_spins();
    auto nbands = mf.get_n_bands();
    auto naos = mf.get_n_aos();

    if ( ov == 'O') // occupied Green's function for negative imaginary time
        tau = - tfg.get_time_nodes()[itau];
    else if ( ov == 'V') // unoccupied Green's function for positive imaginary time
        tau = tfg.get_time_nodes()[itau];
    else
        throw invalid_argument("Invalid occupancy type");
    // alias the real-space Green's function to compute
    auto & gf_Rt = ov == 'O' ? gf_occ_Rt : gf_unocc_Rt;

    // temporary Green's function
    matrix gf_Rt_is_global(naos, naos);

    for (int is = 0; is != nspins; is++)
    {
        gf_Rt_is_global.zero_out();
        auto wg = meanfield.get_weight()[is];
        if ( tau > 0)
            for (int i = 0; i != wg.size; i++)
            {
                wg.c[i] = 1.0 / nkpts * nspins - wg.c[i];
                if (wg.c[i] < 0) wg.c[i] = 0;
            }
        matrix scale(nkpts, nbands);
        // tau-energy phase
        scale = - 0.5 * tau * (mf.get_eigenvals()[is] - mf.get_efermi());
        /* print_matrix("-(e-ef)*tau", scale); */
        for (int ie = 0; ie != scale.size; ie++)
        {
            // NOTE: enforce non-positive phase
            if ( scale.c[ie] > 0) scale.c[ie] = 0;
            scale.c[ie] = std::exp(scale.c[ie]) * wg.c[ie];
        }
        /* print_matrix("exp(-dE*tau)", scale); */
        for (int ik = 0; ik != nkpts; ik++)
        {
            double ang = - klist[ik] * (R * latvec) * TWO_PI;
            complex<double> kphase = complex<double>(cos(ang), sin(ang));
            /* printf("kphase %f %fj\n", kphase.real(), kphase.imag()); */
            auto scaled_wfc_conj = conj(mf.get_eigenvectors()[is][ik]);
            for ( int ib = 0; ib != nbands; ib++)
                LapackConnector::scal(naos, scale(ik, ib), scaled_wfc_conj.c + naos * ib, 1);
            gf_Rt_is_global += (kphase * transpose(mf.get_eigenvectors()[is][ik], false) * scaled_wfc_conj).real();
        }
        if ( tau < 0 ) gf_Rt_is_global *= -1.;
        for ( int I = 0; I != natom; I++)
        {
            const auto I_num = atom_nw[I];
            for (int J = 0; J != natom; J++)
            {
                const size_t J_num = atom_nw[J];
                matrix tmp_green(I_num, J_num);
                for (size_t i = 0; i != I_num; i++)
                {
                    size_t i_glo = atom_iw_loc2glo(I, i);
                    for (size_t j = 0; j != J_num; j++)
                    {
                        size_t j_glo = atom_iw_loc2glo(J, j);
                        tmp_green(i, j) = gf_Rt_is_global(i_glo, j_glo);
                    }
                }
                if (tmp_green.absmax() > gf_R_threshold)
                {
                    // cout<<" max_green_ele:  "<<tmp_green.absmax()<<endl;
                    gf_Rt[is][itau][iR][I][J] = std::move(tmp_green);
                    gf_save++;
                }
                else
                {
                    gf_discard++;
                }
            }
        }
    }
}

void Chi0::build_chi0_q_space_time(const atpair_R_mat_t &LRI_Cs,
                                   const vector<Vector3_Order<int>> &Rlist,
                                   const Vector3_Order<int> &R_period,
                                   const vector<atpair_t> &atpairs_ABF,
                                   const map<int, Vector3_Order<double>> &qlist)
{
    // taus and Rs to compute on MPI task
    // tend to calculate more Rs on one process
    vector<pair<int, int>> itauiRs_local = dispatcher(0, tfg.size(), 0, Rlist.size(), para_mpi.get_myid(), para_mpi.get_size(), true, false);
    map<int, vector<int>> itau2iRs_local;
    for ( auto itauiR: itauiRs_local )
        itau2iRs_local[itauiR.first].push_back(itauiR.second);

    map<int, map<int, atom_mapping<ComplexMatrix>::pair_t_old>> chi0_q_tmp;
    for ( int ifreq = 0; ifreq < tfg.size(); ifreq++ )
        for ( auto iq_qvec: qlist)
        {
            const auto iq = iq_qvec.first;
            for ( auto atpair: atpairs_ABF)
            {
                auto Mu = atpair.first;
                auto Nu = atpair.second;
                chi0_q_tmp[ifreq][iq][Mu][Nu].create(atom_mu[Mu], atom_mu[Nu]);
            }
        }

    omp_lock_t chi0_lock;
    omp_init_lock(&chi0_lock);
    /* if ( atpairs_ABF.size() > itaus_local.size() ) */
    {
        // TODO: add OpenMP paralleling
        for ( auto itau_iRs: itau2iRs_local)
        {
            auto itau = itau_iRs.first;
            auto iRs = itau_iRs.second;
            double tau = tfg.get_freq_nodes()[itau];
            for ( auto atpair: atpairs_ABF)
            {
                atom_t Mu = atpair.first;
                atom_t Nu = atpair.second;
                for ( int is = 0; is != mf.get_n_spins(); is++ )
                {
                    // all Rs is computed at the same time, to avoid repeated calculation of N
                    map<size_t, matrix> chi0_tau;
                    /* if (itau == 0) // debug first itau */
                        chi0_tau = compute_chi0_munu_tau_LRI(gf_occ_Rt[is][itau], gf_unocc_Rt[is][itau],
                                                             LRI_Cs,
                                                             Rlist, R_period, iRs,
                                                             Mu, Nu);
                    /* else continue; // debug first itau */
                    omp_set_lock(&chi0_lock);
                    for ( int iR = 0; iR != Rlist.size(); iR++ )
                    {
                       if ( chi0_tau.count(iR) )
                           for ( int ifreq = 0; ifreq != tfg.size(); ifreq++ )
                           {
                               double trans = tfg.get_costrans_t2f()(ifreq, itau);
                               for ( auto iq_qvec: qlist )
                               {
                                   auto iq = iq_qvec.first;
                                   double arg = iq_qvec.second * (Rlist[iR] * latvec) * TWO_PI;
                                   chi0_q_tmp[ifreq][iq][Mu][Nu] += ComplexMatrix(chi0_tau[iR]) *
                                       (trans * complex<double>(cos(arg), sin(arg)));
                               }
                           }
                    }
                    omp_unset_lock(&chi0_lock);
                }
            }
        }
    }
    /* else */
    /* { */
    /*     throw logic_error("Not implemented"); */
    /* } */

    omp_destroy_lock(&chi0_lock);
    // Reduce to MPI
    for ( int ifreq = 0; ifreq < tfg.size(); ifreq++ )
        for ( auto iq_qvec: qlist )
            for ( auto atpair: atpairs_ABF)
            {
                auto iq = iq_qvec.first;
                auto Mu = atpair.first;
                auto Nu = atpair.second;
                chi0_q[ifreq][iq][Mu][Nu].create(atom_mu[Mu], atom_mu[Nu]);
                /* cout << "nr/nc chi0_q_tmp: " << chi0_q_tmp[ifreq][iq][Mu][Nu].nr << ", "<< chi0_q_tmp[ifreq][iq][Mu][Nu].nc << endl; */
                /* cout << "nr/nc chi0_q: " << chi0_q[ifreq][iq][Mu][Nu].nr << ", "<< chi0_q[ifreq][iq][Mu][Nu].nc << endl; */
                para_mpi.reduce_ComplexMatrix(chi0_q_tmp[ifreq][iq][Mu][Nu], chi0_q[ifreq][iq][Mu][Nu]);
                if (para_mpi.get_myid()==0 && Mu == 0 && Nu == 0 && ifreq == 0 && iq == 0)
                    print_complex_matrix("chi0_q ap 0 0, first freq first q", chi0_q[ifreq][iq][Mu][Nu]);
                chi0_q_tmp[ifreq][iq][Mu].erase(Nu);
            }
}

void Chi0::build_chi0_q_conventional(const atpair_R_mat_t &LRI_Cs,
                                     const vector<Vector3_Order<int>> &Rlist,
                                     const Vector3_Order<int> &R_period,
                                     const vector<atpair_t> &atpair_ABF,
                                     const map<int, Vector3_Order<double>> &qlist)
{
    // TODO: low implementation priority
    throw logic_error("Not implemented");
}

/*!
 * @warning ZMY: Since Green's function are currently stored in a real matrix,
 *          the complex conjugate in the equations are in fact not taken.
 *          But the result seems okay for the test cases, compared to FHI-aims.
 *          Definitely should be resolved in the future.
 */
map<size_t, matrix> compute_chi0_munu_tau_LRI(const map<size_t, atom_mapping<matrix>::pair_t_old> &gf_occ_ab_t,
                                              const map<size_t, atom_mapping<matrix>::pair_t_old> &gf_unocc_ab_t,
                                              const atpair_R_mat_t &LRI_Cs,
                                              const vector<Vector3_Order<int>> &Rlist, const Vector3_Order<int> &R_period,
                                              const vector<int> iRs,
                                              atom_t Mu, atom_t Nu)
{
    map<size_t, matrix> chi0_tau;

    // Store N and N* so one does not need to recompute for each R of chi
    // the extra memory consumption scales as Cs/nR/natoms,
    // which is a small amount compared to Cs for either small (small natoms, large nR) or large (small nR, large natoms) system
    // N(tau) at R and atom K, at(iR, iK)
    map<pair<Vector3_Order<int>, atom_t>, matrix> N;
    // N*(-tau) at R and atom K, at(iR, iK)
    map<pair<Vector3_Order<int>, atom_t>, matrix> Ncnt; // c -> conjugate, nt -> negative time

    const atom_t iI = Mu;
    const atom_t iJ = Nu;
    const size_t n_i = atom_nw[iI];
    const size_t n_mu = atom_mu[Mu];
    const size_t n_j = atom_nw[iJ];
    const size_t n_nu = atom_mu[Nu];

    // pre-compute N and N*
    // only need to calculate N on R-R1, with R in Green's function and R1 from Cs
    // TODO: may reuse the N code to compute N*, by parsing G(t) and G(-t) in same map as done by Shi Rong
    // TODO: decouple N and N*, reduce memory usage
    cout << "Mu: " << Mu << ", Nu: " << Nu << endl;
    for ( auto iR: iRs )
    {
        /* cout << "iR: " << iR << endl; */
        for ( auto const & iK_R1_Cs: LRI_Cs.at(Mu) )
        {
            auto iK = iK_R1_Cs.first;
            /* cout << "iK: " << iK << endl; */
            /* cout << "iR: " << iR << ", iK: " << iK << endl; */
            const size_t n_k = atom_nw[iK];
            for ( auto const & R1_Cs: iK_R1_Cs.second )
            {
                auto const & R1 = R1_Cs.first;
                auto const & RmR1 = Vector3_Order<int>(Rlist[iR] - R1);
                /* cout << "iRmR1: " << iRmR1 << endl; */
                /* cout << "Testing N.count(iRmR1, iK): " << N.count({iRmR1, iK}) << endl; */
                /* cout << "Testing gf_occ_ab_t.count(iR): " << gf_occ_ab_t.count(iR) << endl; */
                // N part, do not recompute and should have the GF counterpart
                if ( N.count({RmR1, iK}) == 0 && gf_occ_ab_t.count(iR) )
                {
                    N[{RmR1, iK}].create(n_nu, n_j*n_k);
                    matrix & N_iRiK = N[{RmR1, iK}];
                    for (auto const & iL_R2_Cs: LRI_Cs.at(Nu))
                    {
                        auto const & iL = iL_R2_Cs.first;
                        const size_t n_l = atom_nw[iL];
                        for ( auto const & R2_Cs: iL_R2_Cs.second )
                        {
                            auto const & R2 = R2_Cs.first;
                            auto mRpR1mR2 = Vector3_Order<int>(-RmR1-R2) % R_period;
                            int i_mRpR1mR2 = get_R_index(Rlist, mRpR1mR2);
                            if ( gf_unocc_ab_t.count(i_mRpR1mR2) )
                            {
                                if (gf_unocc_ab_t.at(i_mRpR1mR2).count(iL) == 0 ||
                                    gf_unocc_ab_t.at(i_mRpR1mR2).at(iL).count(iK) == 0) continue;
                                const matrix & gf_unocc = gf_unocc_ab_t.at(i_mRpR1mR2).at(iL).at(iK);
                                auto const & Cs_nu_jlR2 = R2_Cs.second;
                                assert( Cs_nu_jlR2->nr == n_j * n_l && Cs_nu_jlR2->nc == n_nu );
                                matrix tran_Cs_nu_jlR2 = transpose(*Cs_nu_jlR2);
                                assert( gf_unocc.nr == n_l && gf_unocc.nc == n_k );
                                for ( int i_nu = 0; i_nu != n_nu; i_nu++ )
                                    LapackConnector::gemm('N', 'N', n_j, n_k, n_l, 1.0,
                                                          tran_Cs_nu_jlR2.c+i_nu*n_j*n_l, n_l,
                                                          gf_unocc.c, gf_unocc.nc,
                                                          1.0, N_iRiK.c + i_nu*n_j*n_k, n_k);
                            }
                        }
                    }
                }
                // N* part
                if ( !Ncnt.count({RmR1, iK}) && gf_unocc_ab_t.count(iR) )
                {
                    Ncnt[{RmR1, iK}].create(n_nu, n_j*n_k);
                    matrix & Ncnt_iRiK = Ncnt.at({RmR1, iK});
                    for (auto const & iL_R2_Cs: LRI_Cs.at(Nu))
                    {
                        auto const & iL = iL_R2_Cs.first;
                        const size_t n_l = atom_nw[iL];
                        for ( auto const & R2_Cs: iL_R2_Cs.second )
                        {
                            auto const & R2 = R2_Cs.first;
                            auto mRpR1mR2 = Vector3_Order<int>(-RmR1-R2) % R_period;
                            int i_mRpR1mR2 = get_R_index(Rlist, mRpR1mR2);
                            if ( gf_occ_ab_t.count(i_mRpR1mR2) )
                            {
                                if (gf_occ_ab_t.at(i_mRpR1mR2).count(iL) == 0 ||
                                    gf_occ_ab_t.at(i_mRpR1mR2).at(iL).count(iK) == 0) continue;
                                const matrix & gf_occ = gf_occ_ab_t.at(i_mRpR1mR2).at(iL).at(iK);
                                auto const & Cs_nu_jlR2 = R2_Cs.second;
                                assert( Cs_nu_jlR2->nr == n_j * n_l && Cs_nu_jlR2->nc == n_nu );
                                matrix tran_Cs_nu_jlR2 = transpose(*Cs_nu_jlR2);
                                assert( gf_occ.nr == n_l && gf_occ.nc == n_k );
                                for ( int i_nu = 0; i_nu != n_nu; i_nu++ )
                                    LapackConnector::gemm('N', 'N', n_j, n_k, n_l, 1.0,
                                                          tran_Cs_nu_jlR2.c+i_nu*n_j*n_l, n_l,
                                                          gf_occ.c, n_k,
                                                          1.0, Ncnt_iRiK.c + i_nu*n_j*n_k, n_k);
                            }
                        }
                    }
                }
            }
        }
    }
    cout << "Size N: " << N.size() << endl;
    /* cout << "Precalcualted N, size (iR1 x iK): " << N.size() << endl; */

    for ( auto iR: iRs )
    {
        chi0_tau[iR].create(n_mu, n_nu); // TODO: selectively create by G

        for ( auto const & iK_R1_Cs: LRI_Cs.at(Mu) )
        {
            auto const iK = iK_R1_Cs.first;
            const size_t n_k = atom_nw[iK];
            for (auto const & R1_Cs: iK_R1_Cs.second)
            {
                const auto & Cs = R1_Cs.second;
                // a temporary matrix to store the sum of M(t) and M*(-t)
                matrix MpMc(n_nu, n_i*n_k);
                auto RmR1 = Vector3_Order<int>(Rlist[iR]-R1_Cs.first);
                // M(t)=G(t)N(t)
                if ( gf_occ_ab_t.count(iR) && N.count({RmR1, iK}))
                {
                    if ( gf_occ_ab_t.at(iR).count(iI) != 0 &&
                         gf_occ_ab_t.at(iR).at(iI).count(iJ) != 0)
                    {
                        /* cout << "Computing GN" << endl; */
                        const matrix & gf = gf_occ_ab_t.at(iR).at(iI).at(iJ);
                        const matrix & _N = N.at({RmR1, iK});
                        assert( _N.nr == n_nu && _N.nc == n_j * n_k);
                        for ( int i_nu = 0; i_nu != n_nu; i_nu++ )
                            LapackConnector::gemm('N', 'N', n_i, n_k, n_j, 1.0,
                                                  gf.c, n_j,
                                                  _N.c+i_nu*n_j*n_k, n_k,
                                                  1.0, MpMc.c+i_nu*n_i*n_k, n_k);
                    }
                }
                // M*(-t)=G*(-t)N*(-t)
                if ( gf_unocc_ab_t.count(iR) && N.count({RmR1, iK}))
                {
                    if ( gf_unocc_ab_t.at(iR).count(iI) != 0 &&
                         gf_unocc_ab_t.at(iR).at(iI).count(iJ) != 0)
                    {
                        /* cout << "Computing G*N*" << endl; */
                        const matrix & gf = gf_unocc_ab_t.at(iR).at(iI).at(iJ);
                        const matrix & _Ncnt = Ncnt.at({RmR1, iK});
                        assert( _Ncnt.nr == n_nu && _Ncnt.nc == n_j * n_k);
                        for ( int i_nu = 0; i_nu != n_nu; i_nu++ )
                            LapackConnector::gemm('N', 'N', n_i, n_k, n_j, 1.0,
                                                  gf.c, n_j,
                                                  _Ncnt.c+i_nu*n_j*n_k, n_k,
                                                  1.0, MpMc.c+i_nu*n_i*n_k, n_k);
                    }
                }
                LapackConnector::gemm('T', 'T', n_mu, n_nu, n_i*n_k, 1.0,
                        Cs->c, n_mu,
                        MpMc.c, n_i*n_k, 1.0, chi0_tau[iR].c, n_nu);
            }
        }
    }

    // clean up N to free memory, as they have finished their job
    N.clear();
    Ncnt.clear();
    
    for ( auto iR: iRs )
    {
        matrix X(n_nu, n_j*n_i), Ynt(n_nu, n_i*n_j);
        matrix Xcnt(n_nu, n_j*n_i), Yc(n_nu, n_i*n_j);
        // X
        for ( auto const & iL_R2_Cs: LRI_Cs.at(Nu) )
        {
            auto const iL = iL_R2_Cs.first;
            auto const n_l = atom_nw[iL];
            for ( auto & R2_Cs: iL_R2_Cs.second )
            {
                auto const R2 = R2_Cs.first;
                auto const & Cs = R2_Cs.second;
                auto mR2mR = Vector3_Order<int>(-Rlist[iR]-R2) % R_period;
                auto i_mR2mR = get_R_index(Rlist, mR2mR);
                if ( gf_unocc_ab_t.count(i_mR2mR) &&
                     gf_unocc_ab_t.at(i_mR2mR).count(iL) &&
                     gf_unocc_ab_t.at(i_mR2mR).at(iL).count(iI) )
                {
                    const auto tran_Cs_nu_jl = transpose(*Cs);
                    auto const & gfu = gf_unocc_ab_t.at(i_mR2mR).at(iL).at(iI);
                    for ( int i_nu = 0; i_nu != n_nu; i_nu++ )
                        LapackConnector::gemm('N', 'N', n_j, n_i, n_l, 1.0,
                                              tran_Cs_nu_jl.c+i_nu*n_j*n_l, n_l,
                                              gfu.c, n_i,
                                              1.0, X.c+i_nu*n_j*n_i, n_i);
                }
                if ( gf_occ_ab_t.count(i_mR2mR) &&
                     gf_occ_ab_t.at(i_mR2mR).count(iL) &&
                     gf_occ_ab_t.at(i_mR2mR).at(iL).count(iI) )
                {
                    const auto tran_Cs_nu_jl = transpose(*Cs);
                    // gf shall take conjugate here
                    auto const & gf = gf_occ_ab_t.at(i_mR2mR).at(iL).at(iI);
                    for ( int i_nu = 0; i_nu != n_nu; i_nu++ )
                        LapackConnector::gemm('N', 'N', n_j, n_i, n_l, 1.0,
                                              tran_Cs_nu_jl.c+i_nu*n_j*n_l, n_l,
                                              gf.c, n_i,
                                              1.0, Xcnt.c+i_nu*n_j*n_i, n_i);
                }
            }
        }
        // Y
        for ( auto const & iK_R1_Cs: LRI_Cs.at(Mu) )
        {
            auto const iK= iK_R1_Cs.first;
            auto const n_k = atom_nw[iK];
            for ( auto & R1_Cs: iK_R1_Cs.second )
            {
                auto const R1 = R1_Cs.first;
                auto const Cs = R1_Cs.second;
                auto RmR1 = Vector3_Order<int>(Rlist[iR]-R1) % R_period;
                auto i_RmR1 = get_R_index(Rlist, RmR1);
                if ( gf_occ_ab_t.count(i_RmR1) &&
                     gf_occ_ab_t.at(i_RmR1).count(iK) &&
                     gf_occ_ab_t.at(i_RmR1).at(iK).count(iJ) )
                {
                    const matrix tran_Cs_mu_ik = transpose(*Cs);
                    auto const & gf = gf_occ_ab_t.at(i_RmR1).at(iK).at(iJ);
                    for ( int i_mu = 0; i_mu != n_mu; i_mu++ )
                        // transpose so that result stored in ji order
                        LapackConnector::gemm('T', 'T', n_j, n_i, n_k, 1.0,
                                              gf.c, n_j,
                                              tran_Cs_mu_ik.c+i_mu*n_i*n_k, n_k,
                                              1.0, Ynt.c+i_mu*n_j*n_i, n_i);
                }
                if ( gf_unocc_ab_t.count(i_RmR1) &&
                     gf_unocc_ab_t.at(i_RmR1).count(iK) &&
                     gf_unocc_ab_t.at(i_RmR1).at(iK).count(iJ) )
                {
                    const matrix tran_Cs_mu_ik = transpose(*Cs);
                    auto const & gfu = gf_unocc_ab_t.at(i_RmR1).at(iK).at(iJ);
                    for ( int i_mu = 0; i_mu != n_mu; i_mu++ )
                        // transpose so that result stored in ji order
                        LapackConnector::gemm('T', 'T', n_j, n_i, n_k, 1.0,
                                              gfu.c, n_j,
                                              tran_Cs_mu_ik.c+i_mu*n_i*n_k, n_k,
                                              1.0, Yc.c+i_mu*n_j*n_i, n_i);
                }
            }
        }
        LapackConnector::gemm('N', 'T', n_mu, n_nu, n_i*n_j, 1.0,
                Ynt.c, n_i*n_j,
                X.c, n_i*n_j, 1.0, chi0_tau[iR].c, n_nu);
        LapackConnector::gemm('N', 'T', n_mu, n_nu, n_i*n_j, 1.0,
                Yc.c, n_i*n_j,
                Xcnt.c, n_i*n_j, 1.0, chi0_tau[iR].c, n_nu);
    }
    /* if (para_mpi.get_myid()==0 && Mu == 0 && Nu == 0) */
    if (para_mpi.get_myid()==0 && Mu == 0 && Nu == 1)
    {
        /* for (auto iR_chi0_tau: chi0_tau) */
        /* { */
            /* if (iR_chi0_tau.first!=4) continue; */
            /* cout << Rlist[iR_chi0_tau.first] << " Mu=" << Mu << " Nu=" << Nu; */
            /* print_matrix("chi tauR at first tau and R", iR_chi0_tau.second); */
        /* } */
    }
    /* cout << "Done real-space chi0" << endl; */

    return chi0_tau;
}
