#include <map> 
#include <iostream> 
#include <vector> 
#include "vector3_order.h"
#include "ri.h"
#include "d_chi0.h"
#include "timefreq.h"
#include "chi0.h"
#include "atoms.h"
#include "input.h"
#include "meanfield.h"
#include "chi0.h"
#include "matrix.h"
#include "profiler.h"
 using namespace std;

 vector<Vector3_Order<int>> my_Rlist;
 vector<Vector3_Order<double>> qlist;
 vector<atpair_t> atpairs_ABF;
 Vector3_Order<int> R_period;
 map<int, atom_mapping<map<Vector3_Order<int>, map<double, matrix>>>::pair_t_old> gf_is_R_tau;
 map<int, map<int, double>> rpa_force;

void  d_Chi0::integrate_freq_wc_times_v(map<double, map<Vector3_Order<double>, ComplexMatrix>> &   wc_times_v,
	        	const TFGrids & tfg, map<int, map<Vector3_Order<double>, ComplexMatrix>> & wc_times_v_out) 
  {
    int range_all = N_all_mu;

    for ( auto q_weight: irk_weight)
    {
        qlist.push_back(q_weight.first);
    }


    for (auto it = 0; it != tfg.size(); it++)
    {
    for ( auto q: qlist )
     {
     wc_times_v_out[it][q].create(range_all, range_all);
     wc_times_v_out[it][q].zero_out();      

    for ( int ifreq = 0; ifreq != tfg.size(); ifreq++ )
     {
	auto const freq=tfg.get_freq_nodes()[ifreq];
        const double freq_weight = tfg.find_freq_weight(freq);
        double trans = tfg.get_costrans_t2f()(ifreq, it);
       wc_times_v_out[it][q] +=  freq_weight* trans* wc_times_v[freq][q];     
    }
//     for (auto i=0; i!=range_all; i++) 
//     {	     
//     cout << wc_times_v_out[it][q](i, 1) << endl;	     
//     }
     }	     
     }	    
//  cout << "nanana" << endl;
   
  }; 



   void  d_Chi0::integrate_time_wc_times_v_R_tau(map<int, map<Vector3_Order<double>, ComplexMatrix>> & wc_times_v_out, 
		   const TFGrids & tfg, const MeanField &mf)  
  {
      int range_all = N_all_mu;

      for ( auto R: my_Rlist)
         Rlist_gf.push_back(R);


    vector<pair<int, int>> itauiRs_local = dispatcher(0, tfg.size(), 0, Rlist_gf.size(),
                                                      para_mpi.get_myid(), para_mpi.get_size(), true, false);

     map<int, map<int, ComplexMatrix>>  wc_times_v_R_tau;

    for ( int i=0;i!= itauiRs_local.size();i++)
    {
        auto itau = itauiRs_local[i].first;
        auto tau = tfg.get_time_nodes()[itau];
        auto iR = itauiRs_local[i].second;
        auto R = Rlist_gf[iR];
	wc_times_v_R_tau[iR][itau].create(range_all, range_all);
	wc_times_v_R_tau[iR][itau].zero_out(); 
        for ( auto q: qlist )
	{
          double arg = q * (R * latvec) * TWO_PI;                             
          const complex<double> kphase = complex<double>(cos(arg), sin(arg));
	  wc_times_v_R_tau[iR][itau] += irk_weight[q] * wc_times_v_out[itau][q] * kphase;  
	}
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


     atom_mapping< map<int, map<int, ComplexMatrix>>>::pair_t_old new_wc;

    for ( int i=0;i!= itauiRs_local.size();i++)
    {
        auto itau = itauiRs_local[i].first;
        auto tau = tfg.get_time_nodes()[itau];
        auto iR = itauiRs_local[i].second;
        auto R = Rlist_gf[iR];
        for ( int I = 0; I != natom; I++)
        {
        size_t n_mu=atom_mu[I];
        for ( int J = 0; J != natom; J++)
        {
        size_t n_nu=atom_mu[J];
        new_wc[I][J][iR][itau].create(n_mu, n_nu);   
        new_wc[I][J][iR][itau].zero_out();   
       for (size_t mu = 0; mu != n_mu; ++mu)                                           
        {
       for (size_t nu = 0; nu != n_nu; ++nu)
        {
        new_wc[I][J][iR][itau](mu, nu) = wc_times_v_R_tau[iR][itau](part_range[I] + mu, part_range[J] + nu);
        }                                                                                 
        }

	}

	}

    }


     
  //cout << "bund\n"; 
    //d_Chi0::compute_wc_times_O(wc_times_v_R_tau, itauiRs_local, tfg, mf);
    d_Chi0::compute_wc_times_O(new_wc, itauiRs_local, tfg, mf);
  //cout << "bund\n"; 
  } 


   void d_Chi0::compute_wc_times_O(atom_mapping< map<int, map<int, ComplexMatrix>>>::pair_t_old & wc_times_v_R_tau,
		           vector<pair<int,int>> & itauiRs_local, const TFGrids & tfg, const MeanField &mf)
   {
    

        map<int, map<int, double>> rpa_force_dC;
        ComplexMatrix  O_times_wc;	  
        ComplexMatrix  O_times_wc1;	  
        ComplexMatrix  O_bar_times_wc;	  
        ComplexMatrix  O_bar_times_wc1;	  
        double  constant_1 = 51.422;



        for (auto i_atom=1; i_atom<=natom; i_atom++)
         {
         for (int i_coord=1; i_coord<=3; i_coord++)
         {
         rpa_force_dC[i_coord][i_atom]=0.;
      	 }

	 }




     for ( auto atpair: atpairs_ABF)
     {

            atom_t Mu = atpair.first;
            atom_t Nu = atpair.second;
//    for (auto I_Cs :Cs )
//    { 
//        auto  Mu =I_Cs.first;	    
//    for (auto J_Cs :Cs )
//    {    
//	 auto Nu =J_Cs.first;  
          const size_t i_num = atom_nw[Mu];
          const size_t j_num = atom_nw[Nu];
          const size_t mu_num = atom_mu[Mu];
          const size_t nu_num = atom_mu[Nu];


     for ( auto  K_Cs :Cs[Mu])
     {     
            atom_t Ku = K_Cs.first;
          const size_t k_num = atom_nw[Ku];
     for ( auto  R1_index : K_Cs.second)
      {
          auto   R1 = R1_index.first;
     
       O_times_wc.create(i_num*k_num, mu_num);
       O_times_wc.zero_out();
      if (Mu!=Nu && Nu!=Ku) 
       {
       O_times_wc1.create(j_num*k_num, nu_num);
       O_times_wc1.zero_out();
       }	       

     for ( int i=0;i!= itauiRs_local.size();i++)
      {
        auto itau = itauiRs_local[i].first;
        auto tau = tfg.get_time_nodes()[itau];
        auto iR = itauiRs_local[i].second;
        auto R = Rlist_gf[iR];

    for ( int is = 0; is != mf.get_n_spins(); is++ )
    { 
 //     cout << "nnnnn" << endl;           
        matrix O_matrix;
        if (Mu!=Ku)
	{   	
	O_matrix =2.0 /mf.get_n_spins() * compute_O_s_nu_tau_R_R1(Cs, R_period, is, Mu, Nu, Ku, tau, R, R1);
        O_times_wc += ComplexMatrix(O_matrix)* wc_times_v_R_tau[Nu][Mu][iR][itau];
	}
       if (Mu!=Nu && Nu!=Ku) 
       {
	O_matrix =2.0 /mf.get_n_spins() * compute_O_s_nu_tau_R_R1(Cs, R_period, is, Nu, Mu, Ku, tau, R, R1);
        O_times_wc1 += ComplexMatrix(O_matrix)* wc_times_v_R_tau[Mu][Nu][iR][itau];
       }	       
    }
    }
        for (int i_coord=1; i_coord<=3; i_coord++)
	{
	  if (Mu!=Ku)
          {		  
            matrix d_Cs_tran(transpose(*d_Cs[Mu][Ku][i_coord][R1]));
            ComplexMatrix   prod_matr;
	    prod_matr.create(mu_num, mu_num);
            prod_matr.zero_out();
            prod_matr = ComplexMatrix(d_Cs_tran)* O_times_wc;  
	    complex<double> trace_mat=trace(prod_matr); 
            rpa_force_dC[i_coord][Mu] +=  (-trace_mat.real());	    
            rpa_force_dC[i_coord][Ku] +=   trace_mat.real();	    
          }
       if (Mu!=Nu && Nu!=Ku) 
         {
            matrix d_Cs_tran1(transpose(*d_Cs[Nu][Ku][i_coord][R1]));
            ComplexMatrix   prod_matr1;
            prod_matr1.create(nu_num, nu_num);
            prod_matr1.zero_out();
            prod_matr1 = ComplexMatrix(d_Cs_tran1)* O_times_wc1;  
     	      complex<double> trace_mat1=trace(prod_matr1); 
            rpa_force_dC[i_coord][Nu] +=  (-trace_mat1.real());	    
            rpa_force_dC[i_coord][Ku] +=   trace_mat1.real();	    

	 }

	}

     }



     }

//  Now we take the derivative due to second Cs in Chi.

     for ( auto  L_Cs :Cs[Nu])
     {     
            atom_t Lu = L_Cs.first;
          const size_t l_num = atom_nw[Lu];
     for ( auto  R2_index : L_Cs.second)
      {
          auto   R2 = R2_index.first;

       O_bar_times_wc.create(nu_num, j_num*l_num);
       O_bar_times_wc.zero_out();
      if (Mu!=Nu && Mu!=Lu) 
       {
       O_bar_times_wc1.create(mu_num, i_num*l_num);
       O_bar_times_wc1.zero_out();
       }	       

     for ( int i=0;i!= itauiRs_local.size();i++)
      {
        auto itau = itauiRs_local[i].first;
        auto tau = tfg.get_time_nodes()[itau];
        auto iR = itauiRs_local[i].second;
        auto R = Rlist_gf[iR];

    for ( int is = 0; is != mf.get_n_spins(); is++ )
    { 
       matrix  O_bar_matrix;
      if (Nu!=Lu)
	{   	 
	// O_bar(j*l, mu)	
	O_bar_matrix =2.0 /mf.get_n_spins() * compute_O_bar_s_nu_tau_R_R2(Cs, R_period, is, Mu, Nu, Lu, tau, R, R2);
        O_bar_times_wc += wc_times_v_R_tau[Nu][Mu][iR][itau] *ComplexMatrix(transpose(O_bar_matrix));
       }
       if (Mu!=Nu && Mu!=Lu) 
       {
      // O_bar(i*l, nu)	
   	O_bar_matrix =2.0 /mf.get_n_spins() * compute_O_bar_s_nu_tau_R_R2(Cs, R_period, is, Nu, Mu, Lu, tau, R, R2);
      O_bar_times_wc1 += wc_times_v_R_tau[Mu][Nu][iR][itau] * ComplexMatrix(transpose(O_bar_matrix)) ;
        }	       

    }

      }

        for (int i_coord=1; i_coord<=3; i_coord++)
	{
	  if (Nu!=Lu)
          {		  
            ComplexMatrix   prod_matr2;
            prod_matr2.create(nu_num, nu_num);
            prod_matr2.zero_out();
            prod_matr2 =  O_bar_times_wc * ComplexMatrix(*d_Cs[Nu][Lu][i_coord][R2]);  
	    complex<double> trace_mat2=trace(prod_matr2); 
            rpa_force_dC[i_coord][Nu] +=  (-trace_mat2.real());	    
            rpa_force_dC[i_coord][Lu] +=   trace_mat2.real();	    
          }
       if (Mu!=Nu && Mu!=Lu) 
         {
            ComplexMatrix   prod_matr3;
     	    prod_matr3.create(mu_num, mu_num);
            prod_matr3.zero_out();
            prod_matr3 = O_bar_times_wc1* ComplexMatrix(*d_Cs[Mu][Lu][i_coord][R2]);  
     	    complex<double> trace_mat3=trace(prod_matr3); 
            rpa_force_dC[i_coord][Mu] +=  (-trace_mat3.real());	    
            rpa_force_dC[i_coord][Lu] +=   trace_mat3.real();	    
     	 }
        }


      }
     }     




   }
//   }

         cout << "rpa force dC is " << endl;

        for (auto i_atom=0; i_atom<natom; i_atom++)
	{	
        cout << "atom   "  << i_atom << "  "  <<rpa_force_dC[1][i_atom]*constant_1 << "  " << rpa_force_dC[2][i_atom]*constant_1 << "  " << rpa_force_dC[3][i_atom]*constant_1 << endl; 
	for (int i_coord=1; i_coord<=3; i_coord++)
	{	
	rpa_force[i_coord][i_atom+1] +=(-rpa_force_dC[i_coord][i_atom]);
	}
	}  

	cout << "dV + dC force" << endl;

        for (auto i_atom=1; i_atom<=natom; i_atom++)
	{	
        cout << "atom   "  << i_atom << "  "  <<rpa_force[1][i_atom]*constant_1 << "  " << rpa_force[2][i_atom]*constant_1 << "  " << rpa_force[3][i_atom]*constant_1 << endl; 
	}



 // cout << "bund1\n"; 
   }


matrix d_Chi0::compute_O_s_nu_tau_R_R1(const atpair_R_mat_t &LRI_Cs,
                                       const Vector3_Order<int> &R_period, 
                                       int spin_channel,
                                       atom_t Mu, atom_t Nu, atom_t Ku, double tau, 
				       Vector3_Order<int> R, Vector3_Order<int> R1)
   {

    assert ( tau > 0 );
    // Local RI requires
    const atom_t I_index = Mu;
    const atom_t J_index = Nu;
    const atom_t K_index = Ku;

    const size_t i_num = atom_nw[Mu];
    const size_t mu_num = atom_mu[Mu];

    const size_t j_num = atom_nw[Nu];
    const size_t nu_num = atom_mu[Nu];

    const size_t k_num = atom_nw[Ku];

    int flag_G_IJRt = 0;
    int flag_G_IJRNt = 0;


    const auto & gf_R_tau = gf_is_R_tau.at(spin_channel);
    if (gf_R_tau.at(I_index).count(J_index))
        if (gf_R_tau.at(I_index).at(J_index).count(R))
        {
            if (gf_R_tau.at(I_index).at(J_index).at(R).count(tau))
                flag_G_IJRt = 1;
            if (gf_R_tau.at(I_index).at(J_index).at(R).count(-tau))
                flag_G_IJRNt = 1;
        }


    matrix X_R2(i_num, j_num * nu_num);
    matrix X_conj_R2(i_num, j_num * nu_num);

    for (const auto &L_pair : Cs[J_index])
    {
        const auto L_index = L_pair.first;
        const size_t l_num = atom_nw[L_index];
        if (gf_R_tau.at(I_index).count(L_index))
        {
            for (const auto &R2_index : L_pair.second)
            {
                const auto R2 = R2_index.first;
                const auto &Cs_mat2 = R2_index.second;

                Vector3_Order<int> R_temp_2(Vector3_Order<int>(R2 + R) % R_period);

                if (gf_R_tau.at(I_index).at(L_index).count(R_temp_2))
                {
                    assert(j_num * l_num == (*Cs_mat2).nr);
                    matrix Cs2_reshape(reshape_Cs(j_num, l_num, nu_num, Cs_mat2));
                    
                    if (gf_R_tau.at(I_index).at(L_index).at(R_temp_2).count(tau))
                    {
                        X_R2 += gf_R_tau.at(I_index).at(L_index).at(R_temp_2).at(tau) * Cs2_reshape;
                    }

                    if (gf_R_tau.at(I_index).at(L_index).at(R_temp_2).count(-tau))
                    {
                        X_conj_R2 += gf_R_tau.at(I_index).at(L_index).at(R_temp_2).at(-tau) * Cs2_reshape;
                    }
                    
                }
            }
        }
    }
    matrix X_R2_rs(reshape_mat(i_num, j_num, nu_num, X_R2));
    matrix X_conj_R2_rs(reshape_mat(i_num, j_num, nu_num, X_conj_R2));


                matrix O(i_num, k_num * nu_num);
                matrix Z(k_num, i_num * nu_num);

                if (flag_G_IJRt || flag_G_IJRNt)
                {
                    matrix N_R2(k_num, j_num * nu_num);
                    matrix N_conj_R2(k_num, j_num * nu_num);
                    for (const auto &L_pair : Cs[J_index])
                    {
                        const auto L_index = L_pair.first;
                        const size_t l_num = atom_nw[L_index];
                        if (gf_R_tau.at(K_index).count(L_index))
                        {
                            for (const auto &R2_index : L_pair.second)
                            {
                                const auto R2 = R2_index.first;
                                const auto &Cs_mat2 = R2_index.second;
                                Vector3_Order<int> R_temp_1(Vector3_Order<int>(R + R2 - R1) % R_period);
                                if (gf_R_tau.at(K_index).at(L_index).count(R_temp_1))
                                {
                                    
                                    assert(j_num * l_num == (*Cs_mat2).nr);
                                    matrix Cs2_reshape(reshape_Cs(j_num, l_num, nu_num, Cs_mat2));
                                    if (flag_G_IJRNt && gf_R_tau.at(K_index).at(L_index).at(R_temp_1).count(tau))
                                    {
                                        N_R2 += gf_R_tau.at(K_index).at(L_index).at(R_temp_1).at(tau) * Cs2_reshape;
                                    }
                                    if (flag_G_IJRt && gf_R_tau.at(K_index).at(L_index).at(R_temp_1).count(-tau))
                                    {
                                        N_conj_R2 += gf_R_tau.at(K_index).at(L_index).at(R_temp_1).at(-tau) * Cs2_reshape;
                                    }
                                    
                                }
                            }
                        }
                    }
                    if (flag_G_IJRt)
                    {
                        matrix N_conj_R2_rs(reshape_mat(k_num, j_num, nu_num, N_conj_R2));
                        O += gf_R_tau.at(I_index).at(J_index).at(R).at(tau) * N_conj_R2_rs;
                    }
                    if (flag_G_IJRNt)
                    {
                        matrix N_R2_rs(reshape_mat(k_num, j_num, nu_num, N_R2));
                        O += gf_R_tau.at(I_index).at(J_index).at(R).at(-tau) * N_R2_rs;
                    }
                }
                Vector3_Order<int> R_temp_3(Vector3_Order<int>(R - R1) % R_period);
                if (gf_R_tau.at(K_index).count(J_index))
                {
                    if (gf_R_tau.at(K_index).at(J_index).count(R_temp_3))
                    {
                        if (gf_R_tau.at(K_index).at(J_index).at(R_temp_3).count(-tau))
                        {
                            Z += gf_R_tau.at(K_index).at(J_index).at(R_temp_3).at(-tau) * X_R2_rs;
                        }

                        if (gf_R_tau.at(K_index).at(J_index).at(R_temp_3).count(tau))
                        {
                            Z += gf_R_tau.at(K_index).at(J_index).at(R_temp_3).at(tau) * X_conj_R2_rs;
                        }
                    }
                }

                matrix Z_rs(reshape_mat(k_num, i_num, nu_num, Z));

                O += Z_rs;
                matrix O_matr(reshape_mat_21(i_num, k_num, nu_num, O));


//     cout << "compute_O" << endl;
    return O_matr;
   } 


matrix d_Chi0::compute_O_bar_s_nu_tau_R_R2(const atpair_R_mat_t &LRI_Cs,
                                       const Vector3_Order<int> &R_period, 
                                       int spin_channel,
                                       atom_t Mu, atom_t Nu, atom_t Lu, double tau, 
				       Vector3_Order<int> R, Vector3_Order<int> R2)


{
    prof.start("cal_chi0_element");

    assert ( tau > 0 );
    // Local RI requires
    const atom_t I_index = Mu;
    const atom_t J_index = Nu;
    const atom_t L_index = Lu;

    const size_t i_num = atom_nw[Mu];
    const size_t mu_num = atom_mu[Mu];

    const size_t j_num = atom_nw[Nu];
    const size_t nu_num = atom_mu[Nu];

    const size_t l_num = atom_nw[Lu];


    int flag_G_IJRt = 0;
    int flag_G_IJRNt = 0;
    const auto & gf_R_tau = gf_is_R_tau.at(spin_channel);
    if (gf_R_tau.at(I_index).count(J_index))
        if (gf_R_tau.at(I_index).at(J_index).count(R))
        {
            if (gf_R_tau.at(I_index).at(J_index).at(R).count(tau))
                flag_G_IJRt = 1;

            if (gf_R_tau.at(I_index).at(J_index).at(R).count(-tau))
                flag_G_IJRNt = 1;
        }

    matrix X_R1(j_num, i_num * mu_num);
    matrix X_conj_R1(j_num, i_num * mu_num);

    for (const auto &K_pair : Cs[I_index])
    {
        const auto K_index = K_pair.first;
        const size_t k_num = atom_nw[K_index];
        if (gf_R_tau.at(K_index).count(J_index))
        {
            for (const auto &R1_index : K_pair.second)
            {
                const auto R1 = R1_index.first;
                const auto &Cs_mat1 = R1_index.second;

                Vector3_Order<int> R_temp_2(Vector3_Order<int>(R - R1) % R_period);

                if (gf_R_tau.at(K_index).at(J_index).count(R_temp_2))
                {
                    assert(i_num * k_num == (*Cs_mat1).nr);
                    matrix Cs2_reshape(reshape_Cs(i_num, k_num, mu_num, Cs_mat1));
                    
                    if (gf_R_tau.at(K_index).at(J_index).at(R_temp_2).count(tau))
                    {
                        X_R1 += transpose(gf_R_tau.at(K_index).at(J_index).at(R_temp_2).at(tau)) * Cs2_reshape;
                    }

                    if (gf_R_tau.at(K_index).at(J_index).at(R_temp_2).count(-tau))
                    {
                        X_conj_R1 += transpose(gf_R_tau.at(K_index).at(J_index).at(R_temp_2).at(-tau)) * Cs2_reshape;
                    }
                    
                }
            }
        }
    }
    matrix X_R1_rs(reshape_mat(j_num, i_num, mu_num, X_R1));
    matrix X_conj_R1_rs(reshape_mat(j_num, i_num, mu_num, X_conj_R1));
    
                matrix O(j_num, l_num * mu_num);
                matrix Z(l_num, j_num * mu_num);

                if (flag_G_IJRt || flag_G_IJRNt)
                {
                    matrix N_R1(l_num, i_num * mu_num);
                    matrix N_conj_R1(l_num, i_num * mu_num);
                    for (const auto &K_pair : Cs[I_index])
                    {
                        const auto K_index = K_pair.first;
                        const size_t k_num = atom_nw[K_index];
                        if (gf_R_tau.at(K_index).count(L_index))
                        {
                            for (const auto &R1_index : K_pair.second)
                            {
                                const auto R1 = R1_index.first;
                                const auto &Cs_mat1 = R1_index.second;
                                Vector3_Order<int> R_temp_1(Vector3_Order<int>(R + R2 - R1) % R_period);
                                if (gf_R_tau.at(K_index).at(L_index).count(R_temp_1))
                                {
                                    
                                    assert(i_num * k_num == (*Cs_mat1).nr);
                                    matrix Cs2_reshape(reshape_Cs(i_num, k_num, mu_num, Cs_mat1));
                                    if (flag_G_IJRNt && gf_R_tau.at(K_index).at(L_index).at(R_temp_1).count(tau))
                                    {
                                        N_R1 += transpose(gf_R_tau.at(K_index).at(L_index).at(R_temp_1).at(tau)) * Cs2_reshape;
                                    }
                                    if (flag_G_IJRt && gf_R_tau.at(K_index).at(L_index).at(R_temp_1).count(-tau))
                                    {
                                        N_conj_R1 += transpose(gf_R_tau.at(K_index).at(L_index).at(R_temp_1).at(-tau)) * Cs2_reshape;
                                    }
                                    
                                }
                            }
                        }
                    }
                    if (flag_G_IJRt)
                    {
                        matrix N_conj_R1_rs(reshape_mat(l_num, i_num, mu_num, N_conj_R1));
                        O += transpose(gf_R_tau.at(I_index).at(J_index).at(R).at(tau)) * N_conj_R1_rs;
                    }
                    if (flag_G_IJRNt)
                    {
                        matrix N_R1_rs(reshape_mat(l_num, i_num, mu_num, N_R1));
                        O += transpose(gf_R_tau.at(I_index).at(J_index).at(R).at(-tau)) * N_R1_rs;
                    }
                }
                Vector3_Order<int> R_temp_3(Vector3_Order<int>(R + R2) % R_period);
                if (gf_R_tau.at(I_index).count(L_index))
                {
                    if (gf_R_tau.at(I_index).at(L_index).count(R_temp_3))
                    {
                        if (gf_R_tau.at(I_index).at(L_index).at(R_temp_3).count(-tau))
                        {
                            Z += transpose(gf_R_tau.at(I_index).at(L_index).at(R_temp_3).at(-tau)) * X_R1_rs;
                        }

                        if (gf_R_tau.at(I_index).at(L_index).at(R_temp_3).count(tau))
                        {
                            Z += transpose(gf_R_tau.at(I_index).at(L_index).at(R_temp_3).at(tau)) * X_conj_R1_rs;
                        }
                    }
                }

                matrix Z_rs(reshape_mat(l_num, j_num, mu_num, Z));

                O += Z_rs;
                matrix O_sum(reshape_mat_21(j_num, l_num, mu_num, O));
    
    return O_sum;
}

