#include <iostream>
#include "parallel_mpi.h"
#include "chi0.h"
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
                   const vector<Vector3_Order<double>> &qlist,
                   TFGrids::GRID_TYPES gt, bool use_space_time)
{
    gf_save = gf_discard = 0;
    // reset chi0_q in case the method was called before
    chi0_q.clear();
    for ( int itau = 0; itau < tfg.size(); itau++ )
        for ( int iq = 0; iq < qlist.size(); iq++ )
            for ( auto atpair: atpairs_ABF)
            {
                auto I = atpair.first;
                auto J = atpair.second;
                chi0_q[itau][iq][I][J].create(atom_mu[I], atom_mu[J]);
            }

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
                                   const vector<Vector3_Order<double>> &qlist)
{
    vector<int> itaus_local = dispatcher(0, tfg.size(), para_mpi.get_myid(), para_mpi.get_size(), false);
    omp_lock_t chi0_lock;
    omp_init_lock(&chi0_lock);
    /* if ( atpairs_ABF.size() > itaus_local.size() ) */
    {
        // tau routing
        // TODO: add OpenMP paralleling
        for ( auto itau: itaus_local)
        {
            double tau = tfg.get_freq_nodes()[itau];
            for ( auto atpair: atpairs_ABF)
            {
                atom_t mu = atpair.first;
                atom_t nu = atpair.second;
                for ( int is = 0; is != mf.get_n_spins(); is++ )
                {
                    // all Rs is computed at the same time, to avoid repeated calculation of N
                    vector<matrix> chi0_tau;
                    chi0_tau = compute_chi0_munu_tau_LRI(gf_occ_Rt[is][itau], gf_unocc_Rt[is][itau],
                                                         LRI_Cs,
                                                         Rlist, R_period,
                                                         mu, nu, tau);
                    omp_set_lock(&chi0_lock);
                    for ( int iR = 0; iR != Rlist.size(); iR++ )
                        for ( int ifreq = 0; ifreq != tfg.size(); ifreq++ )
                        {
                            double trans = tfg.get_costrans_t2f()(ifreq, itau);
                            for ( int iq = 0; iq != qlist.size(); iq++ )
                            {
                                double arg = qlist[iq] * (Rlist[iR] * latvec) * TWO_PI;
                                chi0_q[ifreq][iq][mu][nu] += ComplexMatrix(chi0_tau[iR]) *
                                    (trans * std::exp(complex<double>(cos(arg), sin(arg))));
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
}

void Chi0::build_chi0_q_conventional(const atpair_R_mat_t &LRI_Cs,
                                     const vector<Vector3_Order<int>> &Rlist,
                                     const Vector3_Order<int> &R_period,
                                     const vector<atpair_t> &atpair_ABF,
                                     const vector<Vector3_Order<double>> &qlist)
{
    // TODO: low implementation priority
    throw logic_error("Not implemented");
}

vector<matrix> compute_chi0_munu_tau_LRI(const map<size_t, atom_mapping<matrix>::pair_t_old> &gf_occ_ab_t,
                                         const map<size_t, atom_mapping<matrix>::pair_t_old> &gf_unocc_ab_t,
                                         const atpair_R_mat_t &LRI_Cs,
                                         const vector<Vector3_Order<int>> &Rlist, const Vector3_Order<int> &R_period,
                                         atom_t Mu, atom_t Nu, double tau)
{
    vector<matrix> chi0_tau(Rlist.size());
    for (auto &chi0_tauR: chi0_tau)
        chi0_tauR.create(atom_mu[Mu], atom_mu[Nu]);

    // N(-tau) at atom K
    map<atom_t, vector<matrix>> N;
    // N*(-tau) at atom K
    map<atom_t, vector<matrix>> N_cn;

    const size_t n_i = atom_nw[Mu];
    const size_t n_mu = atom_mu[Mu];
    const size_t n_j = atom_nw[Nu];
    const size_t n_nu = atom_mu[Nu];

    // pre-compute N. X can be calculated for each R.
    for ( auto &kR1Cs: LRI_Cs.at(Mu))
    {
    }

    // TODO: real work here
    return chi0_tau;
}
