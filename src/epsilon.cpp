#include <algorithm>
#include "epsilon.h"
#include "input.h"
#include "parallel_mpi.h"
#include "lapack_connector.h"
#include "ri.h"

CorrEnergy compute_RPA_correlation(const Chi0 &chi0, const atpair_k_cplx_mat_t &coulmat)
{
    CorrEnergy corr;
    printf("Begin cal cRPA , pid:  %d\n", para_mpi.get_myid());
    const auto & mf = chi0.mf;
    
    // freq, q
    auto pi_freq_q_Mu_Nu = compute_Pi_q(chi0, coulmat);
    int range_all = 0;

    for (auto &iat : atom_mu)
    {
        range_all += iat.second;
    }
    
    vector<int> part_range;
    part_range.resize(atom_mu.size());
    part_range[0] = 0;
    int count_range = 0;
    for (int I = 0; I != atom_mu.size() - 1; I++)
    {
        count_range += atom_mu[I];
        part_range[I + 1] = count_range;
    }

    cout << "part_range:" << endl;
    for (int I = 0; I != atom_mu.size(); I++)
    {
        cout << part_range[I] << endl;
    }
    cout << "part_range over" << endl;
    // pi_freq_q contains all atoms
    map<double, map<Vector3_Order<double>, ComplexMatrix>> pi_freq_q;
    
    for (const auto &freq_q_MuNupi : pi_freq_q_Mu_Nu)
    {
        const auto freq = freq_q_MuNupi.first;
    
        for (const auto &q_MuNupi : freq_q_MuNupi.second)
        {
            const auto q = q_MuNupi.first;
            const auto MuNupi = q_MuNupi.second;
            pi_freq_q[freq][q].create(range_all, range_all);
    
            ComplexMatrix pi_munu_tmp(range_all, range_all);
            pi_munu_tmp.zero_out();

            for (const auto &Mu_Nupi : MuNupi)
            {
                const auto Mu = Mu_Nupi.first;
                const auto Nupi = Mu_Nupi.second;
                const size_t n_mu = atom_mu[Mu];
                for (const auto &Nu_pi : Nupi)
                {
                    const auto Nu = Nu_pi.first;
                    const auto pimat = Nu_pi.second;
                    const size_t n_nu = atom_mu[Nu];
    
                    for (size_t mu = 0; mu != n_mu; ++mu)
                    {
                        for (size_t nu = 0; nu != n_nu; ++nu)
                        {
                            pi_munu_tmp(part_range[Mu] + mu, part_range[Nu] + nu) += pimat(mu, nu);
                        }
                    }
                }
            }
            if (para_mpi.chi_parallel_type == Parallel_MPI::parallel_type::ATOM_PAIR)
            {
                para_mpi.reduce_ComplexMatrix(pi_munu_tmp, pi_freq_q.at(freq).at(q));
            }
            else
            {
                pi_freq_q.at(freq).at(q) = std::move(pi_munu_tmp);
            }
        }
    }
    if (para_mpi.get_myid() == 0)
    {
        complex<double> tot_RPA_energy(0.0, 0.0);
        map<Vector3_Order<double>, complex<double>> cRPA_q;
        for (const auto &freq_qpi : pi_freq_q)
        {
            const auto freq = freq_qpi.first;
            const double freq_weight = chi0.tfg.find_freq_weight(freq);
            for (const auto &q_pi : freq_qpi.second)
            {
                const auto q = q_pi.first;
                const auto pimat = q_pi.second;
                complex<double> rpa_for_omega_q(0.0, 0.0);
                ComplexMatrix identity(range_all, range_all);
                ComplexMatrix identity_minus_pi(range_all, range_all);
    
                identity.set_as_identity_matrix();
    
                identity_minus_pi = identity - pi_freq_q[freq][q];
                complex<double> det_for_rpa(1.0, 0.0);
                int info_LU = 0;
                int *ipiv = new int[range_all];
                LapackConnector::zgetrf(range_all, range_all, identity_minus_pi, range_all, ipiv, &info_LU);
                for (int ib = 0; ib != range_all; ib++)
                {
                    if (ipiv[ib] != (ib + 1))
                        det_for_rpa = -det_for_rpa * identity_minus_pi(ib, ib);
                    else
                        det_for_rpa = det_for_rpa * identity_minus_pi(ib, ib);
                }
                delete[] ipiv;
    
                complex<double> trace_pi;
                complex<double> ln_det;
                ln_det = std::log(det_for_rpa);
                trace_pi = trace(pi_freq_q.at(freq).at(q));
                cout << "PI trace vector:" << endl;
                cout << endl;
                rpa_for_omega_q = ln_det + trace_pi;
                cout << " ifreq:" << freq << "      rpa_for_omega_k: " << rpa_for_omega_q << "      lnt_det: " << ln_det << "    trace_pi " << trace_pi << endl;
                cRPA_q[q] += rpa_for_omega_q * freq_weight * irk_weight[q] * double(mf.get_n_spins()) / TWO_PI;
                tot_RPA_energy += rpa_for_omega_q * freq_weight * irk_weight[q] * double(mf.get_n_spins()) / TWO_PI;
            }
        }
    
        for (auto &q_crpa : cRPA_q)
        {
            corr.qcontrib[q_crpa.first] = q_crpa.second;
            cout << q_crpa.first << q_crpa.second << endl;
        }
        cout << "gx_num_" << chi0.tfg.size() << "  tot_RPA_energy:  " << setprecision(8) << tot_RPA_energy << endl;
        corr.value = tot_RPA_energy;
    }
    corr.etype = CorrEnergy::type::RPA;
    return corr;
}

CorrEnergy compute_MP2_correlation(const Chi0 &chi0, const atpair_k_cplx_mat_t &coulmat)
{
    CorrEnergy corr;
    corr.etype = CorrEnergy::type::MP2;
    return corr;
}

map<double, map<Vector3_Order<double>, atom_mapping<ComplexMatrix>::pair_t_old>> compute_Pi_q(const Chi0 &chi0, const atpair_k_cplx_mat_t &coulmat)
{
    map<double, map<Vector3_Order<double>, atom_mapping<ComplexMatrix>::pair_t_old>> pi;
    printf("Begin cal_pi_k , pid:  %d\n", para_mpi.get_myid());
    for (auto const & freq_qJQchi0: chi0.get_chi0_q())
    {
        const double freq = freq_qJQchi0.first;
        for (auto &q_JQchi0 : freq_qJQchi0.second)
        {
            Vector3_Order<double> q = q_JQchi0.first;
            for (auto &JQchi0 : q_JQchi0.second )
            {
                const size_t J = JQchi0.first;
                const size_t J_mu = atom_mu[J];
                for (auto &Qchi0 : JQchi0.second)
                {
                    const size_t Q = Qchi0.first;
                    const size_t Q_mu = atom_mu[Q];
                    auto &chi0_mat = Qchi0.second;
                    for (auto &I_p : Vq)
                    {
                        const size_t I = I_p.first;
                        const size_t I_mu = atom_mu[I];
                        pi[freq][q][I][Q].create(I_mu, Q_mu);
                        if (J != Q)
                            pi[freq][q][I][J].create(I_mu, J_mu);
                    }
                }
            }
        }
    }

    for (auto &freq_p : chi0.get_chi0_q())
    {
        const double freq = freq_p.first;
        const auto chi0_freq = freq_p.second;
        for (auto &k_pair : chi0_freq)
        {
            Vector3_Order<double> ik_vec = k_pair.first;
            auto chi0_freq_k = k_pair.second;
            for (auto &J_p : chi0_freq_k)
            {
                const size_t J = J_p.first;
                for (auto &Q_p : J_p.second)
                {
                    const size_t Q = Q_p.first;
                    auto &chi0_mat = Q_p.second;
                    for (auto &I_p : Vq)
                    {
                        const size_t I = I_p.first;
                        // printf("cal_pi  pid: %d , IJQ:  %d  %d  %d\n", para_mpi.get_myid(), I, J, Q);
                        //  cout<<"         pi_IQ: "<<pi_k.at(freq).at(ik_vec).at(I).at(Q)(0,0)<<"   pi_IJ: "<<pi_k.at(freq).at(ik_vec).at(I).at(J)(0,0);
                        if (I <= J)
                        {
                            // cout << "   1"
                            //      << "  Vq: " << (*Vq.at(I).at(J).at(ik_vec))(0, 0) << endl;
                            pi.at(freq).at(ik_vec).at(I).at(Q) += (*Vq.at(I).at(J).at(ik_vec)) * chi0_mat;
                            // if (freq == first_freq)
                            // {
                            //     complex<double> trace_pi;
                            //     trace_pi = trace(pi_k.at(freq).at(ik_vec).at(I).at(Q));
                            //     cout << " IJQ: " << I << " " << J << " " << Q << "  ik_vec: " << ik_vec << "  trace_pi:  " << trace_pi << endl;
                            //     print_complex_matrix("vq:", (*Vq.at(I).at(J).at(ik_vec)));
                            //     print_complex_matrix("chi0:", chi0_mat);
                            //     print_complex_matrix("pi_mat:", pi_k.at(freq).at(ik_vec).at(I).at(Q));
                            // }
                        }
                        else
                        {
                            // cout << "   2"
                            //      << "  Vq: " << transpose(*Vq.at(J).at(I).at(ik_vec), 1)(0, 0) << endl;
                            pi.at(freq).at(ik_vec).at(I).at(Q) += transpose(*Vq.at(J).at(I).at(ik_vec), 1) * chi0_mat;
                        }

                        if (J != Q)
                        {
                            ComplexMatrix chi0_QJ = transpose(chi0_mat, 1);
                            if (I <= Q)
                            {
                                // cout << "   3"
                                //      << "  Vq: " << (*Vq.at(I).at(Q).at(ik_vec))(0, 0) << endl;
                                pi.at(freq).at(ik_vec).at(I).at(J) += (*Vq.at(I).at(Q).at(ik_vec)) * chi0_QJ;
                            }
                            else
                            {
                                // cout << "   4"
                                //      << "  Vq: " << transpose(*Vq.at(J).at(I).at(ik_vec), 1)(0, 0) << endl;
                                pi.at(freq).at(ik_vec).at(I).at(J) += transpose(*Vq.at(Q).at(I).at(ik_vec), 1) * chi0_QJ;
                            }
                        }
                    }
                }
            }
        }
    }
    /* print_complex_matrix(" first_pi_mat:",pi.at(chi0.tfg.get_freq_nodes()[0]).at({0,0,0}).at(0).at(0)); */
    /* print_complex_matrix("  last_pi_mat:",pi.at(chi0.tfg.get_freq_nodes()[0]).at({0,0,0}).at(natom-1).at(natom-1)); */
    return pi;
}


map<double, atpair_k_cplx_mat_t>
compute_Wc_freq_q(const Chi0 &chi0, const atpair_k_cplx_mat_t &coulmat_eps, atpair_k_cplx_mat_t &coulmat_wc)
{
    map<double, atpair_k_cplx_mat_t> Wc_freq_q;
    int range_all = 0;
    for (auto &iat : atom_mu)
    {
        range_all += iat.second;
    }
    const auto part_range = get_part_range();

    // use q-points as the outmost loop, so that square root of Coulomb will not be recalculated at each frequency point
    vector<Vector3_Order<double>> qpts;
    for ( const auto &qMuNuchi: chi0.get_chi0_q().at(chi0.tfg.get_freq_nodes()[0]))
        qpts.push_back(qMuNuchi.first);

    for (const auto &q: qpts)
    {
        int iq = std::distance(klist.begin(), std::find(klist.begin(), klist.end(), q));
        char fn[80];

        ComplexMatrix Vq_all(range_all, range_all);
        for (const auto &Mu_NuqVq: coulmat_eps)
        {
            auto Mu = Mu_NuqVq.first;
            auto n_mu = atom_mu[Mu];
            for ( auto &Nu_qVq: Mu_NuqVq.second )
            {
                auto Nu = Nu_qVq.first;
                if ( 0 == Nu_qVq.second.count(q) ) continue;
                auto n_nu = atom_mu[Nu];
                for ( int i_mu = 0; i_mu != n_mu; i_mu++ )
                    for ( int i_nu = 0; i_nu != n_nu; i_nu++ )
                    {
                        Vq_all(part_range[Mu] + i_mu, part_range[Nu] + i_nu) = (*Nu_qVq.second.at(q))(i_mu, i_nu);
                        Vq_all(part_range[Nu] + i_nu, part_range[Mu] + i_mu) = conj((*Nu_qVq.second.at(q))(i_mu, i_nu));
                    }
            }
        }
        auto sqrtVq_all = power_hemat(Vq_all, 0.5, false, sqrt_coulomb_threshold);
        // sprintf(fn, "sqrtVq_all_q_%d.mtx", iq);
        // print_complex_matrix_mm(sqrtVq_all, fn, 1e-15);

        // truncated (cutoff) Coulomb
        ComplexMatrix Vqcut_all(range_all, range_all);
        for ( auto &Mu_NuqVq: coulmat_wc )
        {
            auto Mu = Mu_NuqVq.first;
            auto n_mu = atom_mu[Mu];
            for ( auto &Nu_qVq: Mu_NuqVq.second )
            {
                auto Nu = Nu_qVq.first;
                if ( 0 == Nu_qVq.second.count(q) ) continue;
                auto n_nu = atom_mu[Nu];
                for ( int i_mu = 0; i_mu != n_mu; i_mu++ )
                    for ( int i_nu = 0; i_nu != n_nu; i_nu++ )
                    {
                        Vqcut_all(part_range[Mu] + i_mu, part_range[Nu] + i_nu) = (*Nu_qVq.second.at(q))(i_mu, i_nu);
                        Vqcut_all(part_range[Nu] + i_nu, part_range[Mu] + i_mu) = conj((*Nu_qVq.second.at(q))(i_mu, i_nu));
                    }
            }
        }
        auto sqrtVqcut_all = power_hemat(Vqcut_all, 0.5, true, sqrt_coulomb_threshold);
        // sprintf(fn, "sqrtVqcut_all_q_%d.mtx", iq);
        // print_complex_matrix_mm(sqrtVqcut_all, fn, 1e-15);
        sprintf(fn, "Vqcut_all_filtered_q_%d.mtx", iq);
        print_complex_matrix_mm(Vqcut_all, fn, 1e-15);
        // save the filtered truncated Coulomb back to the atom mapping object
        for ( auto &Mu_NuqVq: coulmat_wc )
        {
            auto Mu = Mu_NuqVq.first;
            auto n_mu = atom_mu[Mu];
            for ( auto &Nu_qVq: Mu_NuqVq.second )
            {
                auto Nu = Nu_qVq.first;
                if ( 0 == Nu_qVq.second.count(q) ) continue;
                auto n_nu = atom_mu[Nu];
                for ( int i_mu = 0; i_mu != n_mu; i_mu++ )
                    for ( int i_nu = 0; i_nu != n_nu; i_nu++ )
                        (*Nu_qVq.second.at(q))(i_mu, i_nu) = Vqcut_all(part_range[Mu] + i_mu, part_range[Nu] + i_nu);
            }
        }

        ComplexMatrix chi0fq_all(range_all, range_all);
        for (const auto &freq_qMuNuchi: chi0.get_chi0_q())
        {
            auto freq = freq_qMuNuchi.first;
            auto ifreq = chi0.tfg.get_freq_index(freq);
            auto MuNuchi = freq_qMuNuchi.second.at(q);
            for (const auto &Mu_Nuchi: MuNuchi)
            {
                auto Mu = Mu_Nuchi.first;
                auto n_mu = atom_mu[Mu];
                for ( auto &Nu_chi: Mu_Nuchi.second )
                {
                    auto Nu = Nu_chi.first;
                    auto n_nu = atom_mu[Nu];
                    for ( int i_mu = 0; i_mu != n_mu; i_mu++ )
                        for ( int i_nu = 0; i_nu != n_nu; i_nu++ )
                        {
                            chi0fq_all(part_range[Mu] + i_mu, part_range[Nu] + i_nu) = Nu_chi.second(i_mu, i_nu);
                            chi0fq_all(part_range[Nu] + i_nu, part_range[Mu] + i_mu) = conj(Nu_chi.second(i_mu, i_nu));
                        }
                }
            }
            sprintf(fn, "chi0fq_all_q_%d_freq_%d.mtx", iq, ifreq);
            print_complex_matrix_mm(chi0fq_all, fn, 1e-15);

            ComplexMatrix identity(range_all, range_all);
            identity.set_as_identity_matrix();
            auto eps_fq = identity - sqrtVq_all * chi0fq_all * sqrtVq_all;
            // sprintf(fn, "eps_q_%d_freq_%d.mtx", iq, ifreq);
            // print_complex_matrix_mm(eps_fq, fn, 1e-15);

            // invert the epsilon matrix
            power_hemat_onsite(eps_fq, -1);
            auto wc_all = sqrtVqcut_all * (eps_fq - identity) * sqrtVqcut_all;
            // sprintf(fn, "inveps_q_%d_freq_%d.mtx", iq, ifreq);
            // print_complex_matrix_mm(eps_fq, fn, 1e-15);
            sprintf(fn, "wc_q_%d_freq_%d.mtx", iq, ifreq);
            print_complex_matrix_mm(wc_all, fn, 1e-15);

            // save result to the atom mapping object
            for ( auto &Mu_Nuchi: MuNuchi )
            {
                auto Mu = Mu_Nuchi.first;
                auto n_mu = atom_mu[Mu];
                for ( auto &Nu_chi: Mu_Nuchi.second )
                {
                    auto Nu = Nu_chi.first;
                    auto n_nu = atom_mu[Nu];
                    shared_ptr<ComplexMatrix> wc_ptr = make_shared<ComplexMatrix>();
                    wc_ptr->create(n_mu, n_nu);
                    for ( int i_mu = 0; i_mu != n_mu; i_mu++ )
                        for ( int i_nu = 0; i_nu != n_nu; i_nu++ )
                        {
                            (*wc_ptr)(i_mu, i_nu) = wc_all(part_range[Mu] + i_mu, part_range[Nu] + i_nu);
                        }
                    Wc_freq_q[freq][Mu][Nu][q] = wc_ptr;
                }
            }
        }
    }

    return Wc_freq_q;
}

atpair_R_cplx_mat_t
FT_Vq(const atpair_k_cplx_mat_t &coulmat, vector<Vector3_Order<int>> Rlist)
{
    atpair_R_cplx_mat_t VR;
    char fn[80];
    for (auto R: Rlist)
    {
        auto iteR = std::find(Rlist.cbegin(), Rlist.cend(), R);
        auto iR = std::distance(Rlist.cbegin(), iteR);
        for (const auto &Mu_NuqV: coulmat)
        {
            const auto Mu = Mu_NuqV.first;
            const int n_mu = atom_mu[Mu];
            for (const auto &Nu_qV: Mu_NuqV.second)
            {
                const auto Nu = Nu_qV.first;
                const int n_nu = atom_mu[Nu];
                if(!(VR.count(Mu) &&
                     VR.at(Mu).count(Nu) &&
                     VR.at(Mu).at(Nu).count(R)))
                {
                    VR[Mu][Nu][R] = make_shared<ComplexMatrix>();
                    VR[Mu][Nu][R]->create(n_mu, n_nu);
                }
                auto & pV_MuNuR = VR[Mu][Nu][R];
                for (const auto &q_V: Nu_qV.second)
                {
                    auto q = q_V.first;
                    for (auto q_bz: map_irk_ks[q])
                    {
                        /* if (q != klist[1]) continue; // debug zmy, break the result */
                        /* if (q_bz != klist[5]) continue; // debug zmy, break the result */
                        double ang = - q_bz * (R * latvec) * TWO_PI;
                        /* cout << q_bz << ", " << q << endl; */
                        complex<double> kphase = complex<double>(cos(ang), sin(ang));
                        // FIXME: currently support inverse symmetry only
                        if (q_bz == q)
                        {
                            cout << "Direct:  " << q_bz << " => " << q << ", phase = " << kphase << endl;
                            *pV_MuNuR += (*q_V.second) * kphase;
                        }
                        else
                        {
                            cout << "Inverse: " << q_bz << " => " << q << ", phase = " << kphase << endl;
                            *pV_MuNuR += conj(*q_V.second) * kphase;
                        }
                    }
                    // minyez debug: check hermicity of Vq
                    if (iR == 0)
                    {
                        int iq = std::distance(klist.begin(), std::find(klist.begin(), klist.end(), q));
                        sprintf(fn, "Vq_Mu_%zu_Nu_%zu_iq_%d.mtx", Mu, Nu, iq);
                        print_complex_matrix_mm(*q_V.second, fn);
                    }
                    // end minyez debug
                }
                sprintf(fn, "VR_Mu_%zu_Nu_%zu_iR_%zu.mtx", Mu, Nu, iR);
                print_complex_matrix_mm(*pV_MuNuR, fn);
            }
        }
    }
    // myz debug: check the imaginary part of the coulomb matrix
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
    return VR;
}

map<double, atpair_R_cplx_mat_t>
CT_FT_Wc_freq_q(const map<double, atpair_k_cplx_mat_t> &Wc_freq_q,
                const TFGrids &tfg, vector<Vector3_Order<int>> Rlist)
{
    map<double, atpair_R_cplx_mat_t> Wc_tau_R;
    if (!tfg.has_time_grids())
        throw logic_error("TFGrids object does not have time grids");
    const int ngrids = tfg.get_n_grids();
    for (auto R: Rlist)
    {
        for (int itau = 0; itau != ngrids; itau++)
        {
            auto tau = tfg.get_time_nodes()[itau];
            for (int ifreq = 0; ifreq != ngrids; ifreq++)
            {
                auto freq = tfg.get_freq_nodes()[ifreq];
                auto f2t = tfg.get_costrans_f2t()(itau, ifreq);
                if (Wc_freq_q.count(freq))
                {
                    for (const auto &Mu_NuqWc: Wc_freq_q.at(freq))
                    {
                        const auto Mu = Mu_NuqWc.first;
                        const int n_mu = atom_mu[Mu];
                        for (const auto &Nu_qWc: Mu_NuqWc.second)
                        {
                            const auto Nu = Nu_qWc.first;
                            const int n_nu = atom_mu[Nu];
                            if(!(Wc_tau_R.count(tau) &&
                                 Wc_tau_R.at(tau).count(Mu) &&
                                 Wc_tau_R.at(tau).at(Mu).count(Nu) &&
                                 Wc_tau_R.at(tau).at(Mu).at(Nu).count(R)))
                            {
                                Wc_tau_R[tau][Mu][Nu][R] = make_shared<ComplexMatrix>();
                                Wc_tau_R[tau][Mu][Nu][R]->create(n_mu, n_nu);
                            }
                            auto & WtR = Wc_tau_R[tau][Mu][Nu][R];
                            for (const auto &q_Wc: Nu_qWc.second)
                            {
                                auto q = q_Wc.first;
                                for (auto q_bz: map_irk_ks[q])
                                {
                                    double ang = - q_bz * (R * latvec) * TWO_PI;
                                    complex<double> kphase = complex<double>(cos(ang), sin(ang));
                                    complex<double> weight = kphase * f2t;
                                    if (q == q_bz)
                                        *WtR += (*q_Wc.second) * weight;
                                    else
                                        *WtR += conj(*q_Wc.second) * weight;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    // myz debug: check the imaginary part of the matrix
    // NOTE: if G(R) is real, is W(R) real as well?
    for (const auto & tau_MuNuRWc: Wc_tau_R)
    {
        char fn[80];
        auto tau = tau_MuNuRWc.first;
        auto itau = tfg.get_time_index(tau);
        for (const auto & Mu_NuRWc: tau_MuNuRWc.second)
        {
            auto Mu = Mu_NuRWc.first;
            const int n_mu = atom_mu[Mu];
            for (const auto & Nu_RWc: Mu_NuRWc.second)
            {
                auto Nu = Nu_RWc.first;
                const int n_nu = atom_mu[Nu];
                for (const auto & R_Wc: Nu_RWc.second)
                {
                    auto R = R_Wc.first;
                    auto Wc = R_Wc.second;
                    auto iteR = std::find(Rlist.cbegin(), Rlist.cend(), R);
                    auto iR = std::distance(Rlist.cbegin(), iteR);
                    sprintf(fn, "Wc_Mu_%zu_Nu_%zu_iR_%zu_itau_%d.mtx", Mu, Nu, iR, itau);
                    print_complex_matrix_mm(*Wc, fn);
                }
            }
        }
    }
    // end myz debug
    return Wc_tau_R;
}
