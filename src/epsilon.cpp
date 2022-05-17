#include "epsilon.h"
#include "parallel_mpi.h"
#include "lapack_connector.h"

CorrEnergy compute_RPA_correlation(const Chi0 &chi0, const atpair_k_cplx_mat_t &coulmat)
{
    CorrEnergy corr;
    printf("Begin cal cRPA , pid:  %d\n", para_mpi.get_myid());
    const auto & tfg = chi0.tfg;
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
    map<int, map<Vector3_Order<double>, ComplexMatrix>> pi_freq_q;
    
    for (const auto &ifreq_q_MuNupi : pi_freq_q_Mu_Nu)
    {
        const auto ifreq = ifreq_q_MuNupi.first;
    
        for (const auto &q_MuNupi : ifreq_q_MuNupi.second)
        {
            const auto q = q_MuNupi.first;
            const auto MuNupi = q_MuNupi.second;
            pi_freq_q[ifreq][q].create(range_all, range_all);
    
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
            /* if (parall_type == "atom_pair") */
            /* { */
            /*     para_mpi.reduce_ComplexMatrix(tmp_pi_2D, pi_freq_2D.at(freq).at(kvec_c)); */
            /* } */
            /* else */
            {
                pi_freq_q.at(ifreq).at(q) = std::move(pi_munu_tmp);
            }
        }
    }
    if (para_mpi.get_myid() == 0)
    {
        complex<double> tot_RPA_energy(0.0, 0.0);
        map<Vector3_Order<double>, complex<double>> cRPA_q;
        for (const auto &ifreq_qpi : pi_freq_q)
        {
            const auto ifreq = ifreq_qpi.first;
            const double freq_weight = tfg.get_freq_weights()[ifreq];
            for (const auto &q_pi : ifreq_qpi.second)
            {
                const auto q = q_pi.first;
                const auto pimat = q_pi.second;
                complex<double> rpa_for_omega_q(0.0, 0.0);
                ComplexMatrix identity(range_all, range_all);
                ComplexMatrix identity_minus_pi(range_all, range_all);
    
                identity.set_as_identity_matrix();
    
                identity_minus_pi = identity - pi_freq_q[ifreq][q];
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
                trace_pi = trace(pi_freq_q.at(ifreq).at(q));
                cout << "PI trace vector:" << endl;
                cout << endl;
                rpa_for_omega_q = ln_det + trace_pi;
                cout << " ifreq:" << ifreq << "      rpa_for_omega_k: " << rpa_for_omega_q << "      lnt_det: " << ln_det << "    trace_pi " << trace_pi << endl;
                cRPA_q[q] += rpa_for_omega_q * freq_weight * irk_weight[q] * double(mf.get_n_spins()) / TWO_PI;
                tot_RPA_energy += rpa_for_omega_q * freq_weight * irk_weight[q] * double(mf.get_n_spins()) / TWO_PI;
            }
        }
    
        for (auto &q_crpa : cRPA_q)
        {
            corr.qcontrib[q_crpa.first] = q_crpa.second;
            cout << q_crpa.first << q_crpa.second << endl;
        }
        cout << "gx_num_" << tfg.size() << "  tot_RPA_energy:  " << setprecision(8) << tot_RPA_energy << endl;
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

map<int, map<Vector3_Order<double>, atom_mapping<ComplexMatrix>::pair_t_old>> compute_Pi_q(const Chi0 &chi0, const atpair_k_cplx_mat_t &coulmat)
{
    map<int, map<Vector3_Order<double>, atom_mapping<ComplexMatrix>::pair_t_old>> pi;
    printf("Begin cal_pi_k , pid:  %d\n", para_mpi.get_myid());

    return pi;
}
