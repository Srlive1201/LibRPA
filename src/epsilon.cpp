#include <algorithm>
#include "epsilon.h"
#include "input.h"
#include "params.h"
#include "parallel_mpi.h"
#include "lapack_connector.h"
#include "ri.h"
#include "scalapack_connector.h"
#include <omp.h>



complex<double> compute_pi_det_blacs(ComplexMatrix &loc_piT, const int row_nblk, const int col_nblk, const int desc_pi[9], int *ipiv, int &info)
{
    // int range_all = atom_mu_part_range[natom-1]+atom_mu[natom-1];

    // int desc_pi[9];
    // int loc_row, loc_col, info;
    // int row_nblk=1;
    // int col_nblk=1;
    int one=1;
    int range_all= N_all_mu;
    // para_mpi.set_blacs_mat(desc_pi,loc_row,loc_col,range_all,range_all,row_nblk,col_nblk);
    // int *ipiv = new int [loc_row*10];
    // ComplexMatrix loc_piT(loc_col,loc_row);
    
    // for(int i=0;i!=loc_row;i++)
    // {
    //     int global_row = para_mpi.globalIndex(i,row_nblk,para_mpi.nprow,para_mpi.myprow);
    //     int mu;
    //     int I=atom_mu_glo2loc(global_row,mu);
    //     for(int j=0;j!=loc_col;j++)
    //     {
    //         int global_col = para_mpi.globalIndex(j,col_nblk,para_mpi.npcol,para_mpi.mypcol);
    //         int nu;
    //         int J=atom_mu_glo2loc(global_col,nu);
            
    //         if( global_col == global_row)
    //         {
    //             loc_piT(j,i)=complex<double>(1.0,0.0) - pi_freq_q.at(I).at(J)(mu,nu);
    //         }
    //         else
    //         {
    //             loc_piT(j,i)=-1*  pi_freq_q.at(I).at(J)(mu,nu);
    //         }
            
    //     }
    // }
    int DESCPI_T[9];
    // if(out_pi)
    // {
    //     print_complex_real_matrix("first_pi",pi_freq_q.at(0).at(0));
    //     print_complex_real_matrix("first_loc_piT_mat",loc_piT);
    // }
    
    ScalapackConnector::transpose_desc(DESCPI_T,desc_pi);

   // para_mpi.mpi_barrier();
    //printf("   before LU Myid: %d        Available DOS memory = %ld bytes\n",para_mpi.get_myid(), memavail());
    //printf("   before LU myid: %d  range_all: %d,  loc_mat.size: %d\n",para_mpi.get_myid(),range_all,loc_piT.size);
    pzgetrf_(&range_all,&range_all,loc_piT.c,&one,&one,DESCPI_T,ipiv, &info);
    //printf("   after LU myid: %d\n",para_mpi.get_myid());
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
    if(para_mpi.chi_parallel_type == Parallel_MPI::parallel_type::R_TAU)
    {
        compute_RPA_correlation(chi0,coulmat);
        return corr;
    }
    if (para_mpi.get_myid() == 0)
        printf("Calculating EcRPA with BLACS/ScaLAPACK\n");
    // printf("Calculating EcRPA with BLACS, pid:  %d\n", para_mpi.get_myid());
    const auto & mf = chi0.mf;

    //init blacs pi_mat
    int desc_pi[9];
    int loc_row, loc_col, info;
    int row_nblk=1;
    int col_nblk=1;
    int one=1;
    para_mpi.set_blacs_mat(desc_pi,loc_row,loc_col,N_all_mu,N_all_mu,row_nblk,col_nblk);
    int *ipiv = new int [loc_row*10];
    

    map<double, map<Vector3_Order<double>, ComplexMatrix>> pi_freq_q;
    complex<double> tot_RPA_energy(0.0, 0.0);
    map<Vector3_Order<double>, complex<double>> cRPA_q;
    for (const auto &freq_q_MuNuchi0 : chi0.get_chi0_q())
    {
        const auto freq = freq_q_MuNuchi0.first;
        const double freq_weight = chi0.tfg.find_freq_weight(freq);
        for (const auto &q_MuNuchi0 : freq_q_MuNuchi0.second)
        {
            const auto q = q_MuNuchi0.first;
            auto &MuNuchi0 = q_MuNuchi0.second;
            
            ComplexMatrix loc_piT(loc_col,loc_row);
            complex<double> trace_pi(0.0,0.0);
            for(int Mu=0;Mu!=natom;Mu++)
            {
               // printf(" |process %d,  Mu:  %d\n",para_mpi.get_myid(),Mu);
                const size_t n_mu = atom_mu[Mu];
                atom_mapping<ComplexMatrix>::pair_t_old Vq_row=gather_vq_row_q(Mu,coulmat,q);
               // printf("   |process %d,   vq_row.size: %dn",para_mpi.get_myid(),Vq_row[Mu].size());
                ComplexMatrix loc_pi_rowT=compute_Pi_freq_q_row(q,MuNuchi0,Vq_row,Mu);
               // printf("   |process %d,   compute_pi\n",para_mpi.get_myid());
                ComplexMatrix glo_pi_rowT(N_all_mu, n_mu);
                para_mpi.mpi_barrier();
                para_mpi.allreduce_ComplexMatrix(loc_pi_rowT,glo_pi_rowT);
                
                for(int i_mu=0;i_mu!=n_mu;i_mu++)
                    trace_pi+=glo_pi_rowT(atom_mu_part_range[Mu]+i_mu,i_mu);
                //select glo_pi_rowT to pi_blacs 
                for(int i=0;i!=loc_row;i++)
                {
                    int global_row = para_mpi.globalIndex(i,row_nblk,para_mpi.nprow,para_mpi.myprow);
                    int mu_blacs;
                    int I_blacs=atom_mu_glo2loc(global_row,mu_blacs);
                    if(I_blacs== Mu)
                        for(int j=0;j!=loc_col;j++)
                        {
                            int global_col = para_mpi.globalIndex(j,col_nblk,para_mpi.npcol,para_mpi.mypcol);
                            int nu_blacs;
                            int J_blacs=atom_mu_glo2loc(global_col,nu_blacs);

                            if( global_col == global_row)
                            {
                                loc_piT(j,i)=complex<double>(1.0,0.0) - glo_pi_rowT(atom_mu_part_range[J_blacs] + mu_blacs, nu_blacs);
                            }
                            else
                            {
                                loc_piT(j,i)=-1*  glo_pi_rowT(atom_mu_part_range[J_blacs] + mu_blacs, nu_blacs);
                            }

                        }

                }
                
            }
            
            complex<double> ln_det=compute_pi_det_blacs(loc_piT, row_nblk, col_nblk, desc_pi, ipiv, info);
            //para_mpi.mpi_barrier();
            if(para_mpi.get_myid()==0)
            {
                complex<double> rpa_for_omega_q=trace_pi+ln_det;
                //cout << " ifreq:" << freq << "      rpa_for_omega_k: " << rpa_for_omega_q << "      lnt_det: " << ln_det << "    trace_pi " << trace_pi << endl;
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
            //cout << q_crpa.first << q_crpa.second << endl;
        }
        //cout << "gx_num_" << chi0.tfg.size() << "  tot_RPA_energy:  " << setprecision(8)    <<tot_RPA_energy << endl;
    }
    para_mpi.mpi_barrier();
    corr.value = tot_RPA_energy;

    corr.etype = CorrEnergy::type::RPA;
    return corr;
}

CorrEnergy compute_RPA_correlation(const Chi0 &chi0, const atpair_k_cplx_mat_t &coulmat)
{
    CorrEnergy corr;
    if (para_mpi.get_myid() == 0)
        printf("Calculating EcRPA without BLACS/ScaLAPACK\n");
    // printf("Begin cal cRPA , pid:  %d\n", para_mpi.get_myid());
    const auto & mf = chi0.mf;
    
    // freq, q
    map<double, map<Vector3_Order<double>, atom_mapping<ComplexMatrix>::pair_t_old>> pi_freq_q_Mu_Nu;
    if(para_mpi.chi_parallel_type==Parallel_MPI::parallel_type::ATOM_PAIR)
        pi_freq_q_Mu_Nu = compute_Pi_q_MPI(chi0, coulmat);
    else 
        pi_freq_q_Mu_Nu = compute_Pi_q(chi0, coulmat);

    int range_all = N_all_mu;

    
    vector<int> part_range;
    part_range.resize(atom_mu.size());
    part_range[0] = 0;
    int count_range = 0;
    for (int I = 0; I != atom_mu.size() - 1; I++)
    {
        count_range += atom_mu[I];
        part_range[I + 1] = count_range;
    }

    // cout << "part_range:" << endl;
    // for (int I = 0; I != atom_mu.size(); I++)
    // {
    //     cout << part_range[I] << endl;
    // }
    // cout << "part_range over" << endl;

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
                // cout << "PI trace vector:" << endl;
                // cout << endl;
                rpa_for_omega_q = ln_det + trace_pi;
                // cout << " ifreq:" << freq << "      rpa_for_omega_k: " << rpa_for_omega_q << "      lnt_det: " << ln_det << "    trace_pi " << trace_pi << endl;
                cRPA_q[q] += rpa_for_omega_q * freq_weight * irk_weight[q] * double(mf.get_n_spins()) / TWO_PI;
                tot_RPA_energy += rpa_for_omega_q * freq_weight * irk_weight[q] * double(mf.get_n_spins()) / TWO_PI;
            }
        }
    
        for (auto &q_crpa : cRPA_q)
        {
            corr.qcontrib[q_crpa.first] = q_crpa.second;
        }
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
    // printf("Begin cal_pi_k , pid:  %d\n", para_mpi.get_myid());
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
                    // auto &chi0_mat = Qchi0.second;
                    for (int I=0;I!=natom;I++)
                    {
                        //const size_t I = I_p.first;
                        const size_t I_mu = atom_mu[I];
                        pi[freq][q][I][Q].create(I_mu, Q_mu);
                        if (J != Q)
                            pi[freq][q][I][J].create(I_mu, J_mu);
                    }
                }
            }
            // if(freq==chi0.tfg.get_freq_nodes()[0])
            //     for(auto &Ip:pi[freq][q])
            //         for(auto &Jp:Ip.second)
            //             printf("  |process  %d, pi atpair: %d, %d \n",para_mpi.get_myid(),Ip.first,Jp.first);
        }
        
    }

    // ofstream fp;
    // std::stringstream ss;
    // ss<<"out_pi_rank_"<<para_mpi.get_myid()<<".txt";
    // fp.open(ss.str());
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
                    for (int I=0;I!=natom;I++)
                    {
                        //const size_t I = I_p.first;
                        //printf("cal_pi  pid: %d , IJQ:  %d  %d  %d\n", para_mpi.get_myid(), I, J, Q);
                        //  cout<<"         pi_IQ: "<<pi_k.at(freq).at(ik_vec).at(I).at(Q)(0,0)<<"   pi_IJ: "<<pi_k.at(freq).at(ik_vec).at(I).at(J)(0,0);
                        if (I <= J)
                        {
                            // if (freq == chi0.tfg.get_freq_nodes()[0])
                            //     printf("cal_pi  pid: %d , IJQ:  %d  %d  %d   type: %d \n", para_mpi.get_myid(), I, J, Q,1);
                            //      << "  Vq: " << (*Vq.at(I).at(J).at(ik_vec))(0, 0) << endl;
                            pi.at(freq).at(ik_vec).at(I).at(Q) += (*Vq.at(I).at(J).at(ik_vec)) * chi0_mat;
                            // if (freq == chi0.tfg.get_freq_nodes()[0])
                            // {
                            //     std:stringstream sm;
                            //     complex<double> trace_pi;
                            //     trace_pi = trace(pi.at(freq).at(ik_vec).at(I).at(Q));
                            //     sm << " IJQ: " << I << " " << J << " " << Q << "  ik_vec: " << ik_vec << "  trace_pi:  " << trace_pi << endl;
                            //     print_complex_matrix_file(sm.str().c_str(), (*Vq.at(I).at(J).at(ik_vec)),fp,false);
                            //     print_complex_matrix_file("chi0:", chi0_mat,fp,false);
                            //     print_complex_matrix_file("pi_mat:", pi.at(freq).at(ik_vec).at(I).at(Q),fp,false);
                            // }
                        }
                        else
                        {
                            // if (freq == chi0.tfg.get_freq_nodes()[0])
                            //     printf("cal_pi  pid: %d , IJQ:  %d  %d  %d   type: %d \n", para_mpi.get_myid(), I, J, Q,2);
                            //      << "  Vq: " << transpose(*Vq.at(J).at(I).at(ik_vec), 1)(0, 0) << endl;
                            pi.at(freq).at(ik_vec).at(I).at(Q) += transpose(*Vq.at(J).at(I).at(ik_vec), 1) * chi0_mat;
                        }

                        if (J != Q)
                        {
                            ComplexMatrix chi0_QJ = transpose(chi0_mat, 1);
                            if (I <= Q)
                            {
                                // if (freq == chi0.tfg.get_freq_nodes()[0])
                                //     printf("cal_pi  pid: %d , IJQ:  %d  %d  %d   type: %d \n", para_mpi.get_myid(), I, J, Q,3);
                                //      << "  Vq: " << (*Vq.at(I).at(Q).at(ik_vec))(0, 0) << endl;
                                pi.at(freq).at(ik_vec).at(I).at(J) += (*Vq.at(I).at(Q).at(ik_vec)) * chi0_QJ;
                            }
                            else
                            {
                                // if (freq == chi0.tfg.get_freq_nodes()[0])
                                //     printf("cal_pi  pid: %d , IJQ:  %d  %d  %d   type: %d \n", para_mpi.get_myid(), I, J, Q,4);
                                //      << "  Vq: " << transpose(*Vq.at(J).at(I).at(ik_vec), 1)(0, 0) << endl;
                                pi.at(freq).at(ik_vec).at(I).at(J) += transpose(*Vq.at(Q).at(I).at(ik_vec), 1) * chi0_QJ;
                            }
                        }
                    }
                }
            }
        }
    }
    //fp.close();
    // print_complex_matrix(" first_pi_mat:",pi.at(chi0.tfg.get_freq_nodes()[0]).at({0,0,0}).at(0).at(0)); 
    /* print_complex_matrix("  last_pi_mat:",pi.at(chi0.tfg.get_freq_nodes()[0]).at({0,0,0}).at(natom-1).at(natom-1)); */
    return pi;
}

map<double, map<Vector3_Order<double>, atom_mapping<ComplexMatrix>::pair_t_old>> compute_Pi_q_MPI(const Chi0 &chi0, const atpair_k_cplx_mat_t &coulmat)
{
    map<double, map<Vector3_Order<double>, atom_mapping<ComplexMatrix>::pair_t_old>> pi;
    // printf("Begin cal_pi_k , pid:  %d\n", para_mpi.get_myid());
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
                    // auto &chi0_mat = Qchi0.second;
                    for (int I=0;I!=natom;I++)
                    {
                        //const size_t I = I_p.first;
                        const size_t I_mu = atom_mu[I];
                        pi[freq][q][I][Q].create(I_mu, Q_mu);
                        if (J != Q)
                            pi[freq][q][I][J].create(I_mu, J_mu);
                    }
                }
            }
            // if(freq==chi0.tfg.get_freq_nodes()[0])
            //     for(auto &Ip:pi[freq][q])
            //         for(auto &Jp:Ip.second)
            //             printf("  |process  %d, pi atpair: %d, %d \n",para_mpi.get_myid(),Ip.first,Jp.first);
        }
        
    }

    // ofstream fp;
    // std::stringstream ss;
    // ss<<"out_pi_rank_"<<para_mpi.get_myid()<<".txt";
    // fp.open(ss.str());
    for (auto &k_pair : irk_weight)
    {
        Vector3_Order<double> ik_vec = k_pair.first;
        for (int I=0;I!=natom;I++)
        {
            atom_mapping<ComplexMatrix>::pair_t_old Vq_row=gather_vq_row_q(I,coulmat,ik_vec);
            for (auto &freq_p : chi0.get_chi0_q())
            {
                const double freq = freq_p.first;
                const auto chi0_freq = freq_p.second;
            
                
                auto chi0_freq_k = freq_p.second.at(ik_vec);

                for (auto &J_p : chi0_freq_k)
                {
                    const size_t J = J_p.first;
                    for (auto &Q_p : J_p.second)
                    {
                        const size_t Q = Q_p.first;
                        auto &chi0_mat = Q_p.second;
                    
                        //const size_t I = I_p.first;
                        //printf("cal_pi  pid: %d , IJQ:  %d  %d  %d\n", para_mpi.get_myid(), I, J, Q);
                        //  cout<<"         pi_IQ: "<<pi_k.at(freq).at(ik_vec).at(I).at(Q)(0,0)<<"   pi_IJ: "<<pi_k.at(freq).at(ik_vec).at(I).at(J)(0,0);
                        
                        // if (freq == chi0.tfg.get_freq_nodes()[0])
                        //     printf("cal_pi  pid: %d , IJQ:  %d  %d  %d   type: %d \n", para_mpi.get_myid(), I, J, Q,1);
                        //      << "  Vq: " << (*Vq.at(I).at(J).at(ik_vec))(0, 0) << endl;
                        pi.at(freq).at(ik_vec).at(I).at(Q) += Vq_row.at(I).at(J) * chi0_mat;
                        // if (freq == chi0.tfg.get_freq_nodes()[0])
                        // {
                        //     std:stringstream sm;
                        //     complex<double> trace_pi;
                        //     trace_pi = trace(pi.at(freq).at(ik_vec).at(I).at(Q));
                        //     sm << " IJQ: " << I << " " << J << " " << Q << "  ik_vec: " << ik_vec << "  trace_pi:  " << trace_pi << endl;
                        //     print_complex_matrix_file(sm.str().c_str(), Vq_row.at(I).at(J),fp,false);
                        //     print_complex_matrix_file("chi0:", chi0_mat,fp,false);
                        //     print_complex_matrix_file("pi_mat:", pi.at(freq).at(ik_vec).at(I).at(Q),fp,false);
                        // }  

                        if (J != Q)
                        {
                            ComplexMatrix chi0_QJ = transpose(chi0_mat, 1);
                            // if (freq == chi0.tfg.get_freq_nodes()[0])
                            //     printf("cal_pi  pid: %d , IJQ:  %d  %d  %d   type: %d \n", para_mpi.get_myid(), I, J,Q,3);
                            //      << "  Vq: " << (*Vq.at(I).at(Q).at(ik_vec))(0, 0) << endl;
                            pi.at(freq).at(ik_vec).at(I).at(J) += Vq_row.at(I).at(Q) * chi0_QJ;
                        }
                    }
                }
            }
        }
    }
    //fp.close();
    // print_complex_matrix(" first_pi_mat:",pi.at(chi0.tfg.get_freq_nodes()[0]).at({0,0,0}).at(0).at(0)); 
    /* print_complex_matrix("  last_pi_mat:",pi.at(chi0.tfg.get_freq_nodes()[0]).at({0,0,0}).at(natom-1).at(natom-1)); */
    return pi;
}



ComplexMatrix compute_Pi_freq_q_row(const Vector3_Order<double> &ik_vec, const atom_mapping<ComplexMatrix>::pair_t_old &chi0_freq_q, const atom_mapping<ComplexMatrix>::pair_t_old &Vq_row, const int &I)
{
    map<size_t,ComplexMatrix> pi;
    // printf("Begin cal_pi_k , pid:  %d\n", para_mpi.get_myid());
    auto I_mu=atom_mu[I];
    for(int J=0;J!=natom;J++)
        pi[J].create(I_mu,atom_mu[J]);
   
    for (auto &J_p : chi0_freq_q)
    {
        const size_t J = J_p.first;
        for (auto &Q_p : J_p.second)
        {
            const size_t Q = Q_p.first;
            auto &chi0_mat = Q_p.second;
            pi.at(Q) += Vq_row.at(I).at(J) * chi0_mat;
            if (J != Q)
            {
                ComplexMatrix chi0_QJ = transpose(chi0_mat, 1);
                pi.at(J) += Vq_row.at(I).at(Q) * chi0_QJ;
            }
        }
    }
    ComplexMatrix pi_row(N_all_mu,atom_mu[I]);
    complex<double> *pi_row_ptr=pi_row.c;
    for(auto &Jp:pi)
    {
        auto J=Jp.first;
        auto J_mu=atom_mu[J];
        const auto length=sizeof(complex<double>)* I_mu *J_mu;
        memcpy(pi_row_ptr, pi.at(J).c,length);
        pi_row_ptr+=I_mu *J_mu;
    }
    return pi_row;
}

atom_mapping<ComplexMatrix>::pair_t_old gather_vq_row_q(const int &I, const atpair_k_cplx_mat_t &coulmat, const Vector3_Order<double> &ik_vec)
{
    auto I_mu=atom_mu[I];
    atom_mapping<ComplexMatrix>::pair_t_old Vq_row;
    for(int J_tmp=0;J_tmp!=natom;J_tmp++)
    {
        auto J_mu=atom_mu[J_tmp];
        ComplexMatrix loc_vq(atom_mu[I],atom_mu[J_tmp]);
        Vq_row[I][J_tmp].create(atom_mu[I],atom_mu[J_tmp]);
       // const auto length=sizeof(complex<double>)* I_mu *J_mu;
       // complex<double> *loc_vq_ptr=loc_vq.c;
        if(I<=J_tmp)
        {
            if(Vq.count(I))
                if(Vq.at(I).count(J_tmp))
                    loc_vq=*Vq.at(I).at(J_tmp).at(ik_vec);
        }
        else
        {
            if(Vq.count(J_tmp))
                if(Vq.at(J_tmp).count(I))
                    loc_vq=transpose(*Vq.at(J_tmp).at(I).at(ik_vec), 1);
                    
        }
        para_mpi.allreduce_ComplexMatrix(loc_vq,Vq_row[I][J_tmp]);
    }
    return Vq_row;
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
        auto sqrtVq_all = power_hemat(Vq_all, 0.5, false, params.sqrt_coulomb_threshold);
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
        auto sqrtVqcut_all = power_hemat(Vqcut_all, 0.5, true, params.sqrt_coulomb_threshold);
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
            // const int n_mu = atom_mu[Mu];
            for (const auto & Nu_RWc: Mu_NuRWc.second)
            {
                auto Nu = Nu_RWc.first;
                // const int n_nu = atom_mu[Nu];
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
