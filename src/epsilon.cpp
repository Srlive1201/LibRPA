#include "epsilon.h"
#include "parallel_mpi.h"
#include "lapack_connector.h"
#include "scalapack_connector.h"
complex<double> compute_pi_det(map<size_t, map<size_t, ComplexMatrix>> &pi_freq_q,bool out_pi)
{
    int range_all = atom_mu_part_range[natom-1]+atom_mu[natom-1];

    int desc_pi[9];
    int loc_row, loc_col, info;
    int row_nblk=1;
    int col_nblk=1;
    int one=1;
    para_mpi.set_blacs_mat(desc_pi,loc_row,loc_col,range_all,range_all,row_nblk,col_nblk);
    int *ipiv = new int [loc_row];
    ComplexMatrix loc_piT(loc_col,loc_row);
    
    for(int i=0;i!=loc_row;i++)
    {
        int global_row = para_mpi.globalIndex(i,row_nblk,para_mpi.nprow,para_mpi.myprow);
        int mu;
        int I=atom_mu_glo2loc(global_row,mu);
        for(int j=0;j!=loc_col;j++)
        {
            int global_col = para_mpi.globalIndex(j,col_nblk,para_mpi.npcol,para_mpi.mypcol);
            int nu;
            int J=atom_mu_glo2loc(global_col,nu);
            // if(out_pi)
            // {
            //     cout<<"i,j: "<<i<<" "<<j<<"  glo: ("<<global_row<<", "<<global_col<<" )"<<"    I J mu nu: "<<I<<"  "<<J<<"  "<<mu<<"  "<<nu<<endl;
            // }
            
            if( global_col == global_row)
            {
                loc_piT(j,i)=complex<double>(1.0,0.0) - pi_freq_q.at(I).at(J)(mu,nu);
            }
            else
            {
                loc_piT(j,i)=-1*  pi_freq_q.at(I).at(J)(mu,nu);
            }
            
        }
    }
    int DESCPI_T[9];
    // if(out_pi)
    // {
    //     print_complex_real_matrix("first_pi",pi_freq_q.at(0).at(0));
    //     print_complex_real_matrix("first_loc_piT_mat",loc_piT);
    // }
    ScalapackConnector::transpose_desc(DESCPI_T,desc_pi);
    pzgetrf_(&range_all,&range_all,loc_piT.c,&one,&one,DESCPI_T,ipiv, &info);
    complex<double> ln_det_loc(0.0,0.0);
    complex<double> ln_det_all(0.0,0.0);
    for(int ig=0;ig!=range_all;ig++)
    {
        int locr=para_mpi.localIndex(ig,row_nblk,para_mpi.nprow,para_mpi.myprow);
		int locc=para_mpi.localIndex(ig,col_nblk,para_mpi.npcol,para_mpi.mypcol);
		if(locr >=0 && locc>=0)
		{
			// if(ipiv[locr]!=(ig+1))
			// 	det_loc=-1*det_loc * loc_piT(locc,locr);
			// else
			// 	det_loc=det_loc * loc_piT(locc,locr);
            if(loc_piT(locc,locr).real()>0)
                ln_det_loc+=std::log(loc_piT(locc,locr));
            else
                ln_det_loc+=std::log(-1 * loc_piT(locc,locr));
		}
    }
    MPI_Allreduce(&ln_det_loc,&ln_det_all,1,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD);
    return ln_det_all;
}


CorrEnergy compute_RPA_correlation_blacs(const Chi0 &chi0, const atpair_k_cplx_mat_t &coulmat)
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
    complex<double> tot_RPA_energy(0.0, 0.0);
    map<Vector3_Order<double>, complex<double>> cRPA_q;
    for (const auto &freq_q_MuNupi : pi_freq_q_Mu_Nu)
    {
        const auto freq = freq_q_MuNupi.first;
        const double freq_weight = chi0.tfg.find_freq_weight(freq);
        for (const auto &q_MuNupi : freq_q_MuNupi.second)
        {
            const auto q = q_MuNupi.first;
            auto MuNupi = q_MuNupi.second;
            pi_freq_q[freq][q].create(range_all, range_all);
    
            ComplexMatrix pi_munu_tmp(range_all, range_all);
            pi_munu_tmp.zero_out();

            if (para_mpi.chi_parallel_type == Parallel_MPI::parallel_type::ATOM_PAIR)
            {
                for (int Mu=0;Mu!=atom_mu.size();Mu++)
                {
                    const size_t n_mu = atom_mu[Mu];
                    for (int Nu=0;Nu!=atom_mu.size();Nu++)
                    {
                        const size_t n_nu = atom_mu[Nu];

                        ComplexMatrix loc_pi_Mu_Nu(n_mu,n_nu);
                        ComplexMatrix rd_pi_Mu_Nu(n_mu,n_nu);
                        if(MuNupi.count(Mu))
                        {
                            if(MuNupi.at(Mu).count(Nu))
                                loc_pi_Mu_Nu=MuNupi.at(Mu).at(Nu);
                        }

                        para_mpi.allreduce_ComplexMatrix(loc_pi_Mu_Nu,rd_pi_Mu_Nu);

                        MuNupi[Mu][Nu]=std::move(rd_pi_Mu_Nu);

                    }
                }
            }
            para_mpi.mpi_barrier();
            complex<double> trace_pi(0.0,0.0);
            for (int Mu=0;Mu!=atom_mu.size();Mu++)
            {
                const size_t n_mu = atom_mu[Mu];
                for(int mu=0;mu!=n_mu;mu++)
                    trace_pi+=MuNupi.at(Mu).at(Mu)(mu,mu);
            }
            bool out_pi=false;
            if(freq==chi0.tfg.get_freq_nodes()[0] && q == Vector3_Order<double>({0,0,0}))
                out_pi=true;
            complex<double> ln_det=compute_pi_det(MuNupi,out_pi);
            //para_mpi.mpi_barrier();
            if(para_mpi.get_myid()==0)
            {
                complex<double> rpa_for_omega_q=trace_pi+ln_det;
                cout << " ifreq:" << freq << "      rpa_for_omega_k: " << rpa_for_omega_q << "      lnt_det: " << ln_det << "    trace_pi " << trace_pi << endl;
                cRPA_q[q] += rpa_for_omega_q * freq_weight * irk_weight[q] * double(mf.get_n_spins()) / TWO_PI;
                tot_RPA_energy += rpa_for_omega_q * freq_weight * irk_weight[q] * double(mf.get_n_spins()) / TWO_PI;
            }
        }
    }
    
    if(para_mpi.get_myid()==0)
    {
        for (auto &q_crpa : cRPA_q)
        {
            corr.qcontrib[q_crpa.first] = q_crpa.second;
            cout << q_crpa.first << q_crpa.second << endl;
        }
        cout << "gx_num_" << chi0.tfg.size() << "  tot_RPA_energy:  " << setprecision(8)    <<tot_RPA_energy << endl;
    }
    para_mpi.mpi_barrier();
    corr.value = tot_RPA_energy;

    corr.etype = CorrEnergy::type::RPA;
    return corr;
}


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
    // print_complex_matrix(" first_pi_mat:",pi.at(chi0.tfg.get_freq_nodes()[0]).at({0,0,0}).at(0).at(0)); 
    /* print_complex_matrix("  last_pi_mat:",pi.at(chi0.tfg.get_freq_nodes()[0]).at({0,0,0}).at(natom-1).at(natom-1)); */
    return pi;
}
