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
#include "lapack_connector.h"
#include "scalapack_connector.h"
//#include "parallel_mpi.h"
#include <omp.h>
 using namespace std;

 vector<Vector3_Order<int>> my_Rlist;
 vector<Vector3_Order<double>> qlist;
 vector<atpair_t> atpairs_ABF;
 Vector3_Order<int> R_period;
 map<int, atom_mapping<map<Vector3_Order<int>, map<double, matrix>>>::pair_t_old> gf_is_R_tau;
 map<int, atom_mapping<map<Vector3_Order<int>, map<double, matrix>>>::pair_t_old> d_gf_R_tau;
 map<int, map<int, atom_mapping<map<Vector3_Order<int>, map<double, matrix>>>::pair_t_old>> d_gf_R_tau_full;
 map<int, map<int, double>> rpa_force;
 matrix n_homo_k, n_lumo_k;

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

      cout << "start dG part now" << endl;
      find_n_homo_k_n_lumo_k(mf);


        for ( auto tau: tfg.get_time_nodes() )
        {
            for ( auto R: Rlist_gf)
            {

              build_d_gf_full_unocc_R_plus_tau( R,  tau, mf);
              build_d_gf_full_occ_R_minus_tau ( R, -tau, mf);
            }
        }

// **********************************************************************************************************


        map<int, map<int, double>> rpa_force_dG;

       for (int i_coord=1; i_coord<=3; i_coord++)
       {
       for (int i_atom =0; i_atom !=natom; i_atom ++ )
       {		   
         rpa_force_dG[i_coord][i_atom]=0.;
       }
       }


     for ( auto atpair: atpairs_ABF)
     {

            atom_t Mu = atpair.first;
            atom_t Nu = atpair.second;
          const size_t i_num = atom_nw[Mu];
          const size_t j_num = atom_nw[Nu];
          const size_t mu_num = atom_mu[Mu];
          const size_t nu_num = atom_mu[Nu];

         ComplexMatrix self_energy_times_dG2 (j_num, j_num); 
         ComplexMatrix self_energy_times_dG3 (i_num, i_num); 
     for ( int i=0;i!= itauiRs_local.size();i++)
      {
        auto itau = itauiRs_local[i].first;
        auto tau = tfg.get_time_nodes()[itau];
        auto iR = itauiRs_local[i].second;
        auto R = Rlist_gf[iR];
    for ( int is = 0; is != mf.get_n_spins(); is++ )
    { 
      const auto & gf_R_tau = gf_is_R_tau.at(is);

   map<int, map<int, double>>  rpa_force_dG_kj;  
   map<int, map<int, double>>  rpa_force_dG_jk;  
   map<int, map<int, double>> rpa_force_dG_ji;
   map<int, map<int, double>> rpa_force_dG_ij;
   map<int, map<int, double>> rpa_force_dG_lk;
   map<int, map<int, double>> rpa_force_dG_kl;
   map<int, map<int, double>> rpa_force_dG_il;
   map<int, map<int, double>> rpa_force_dG_li;


   for (int i_coord=1; i_coord<=3; i_coord++)
   {
   for (int i_atom =0; i_atom !=natom; i_atom ++ )
   {		   
    	  rpa_force_dG_kj[i_coord][i_atom] =0.;
    	  rpa_force_dG_jk[i_coord][i_atom] =0.;
          rpa_force_dG_ji[i_coord][i_atom] =0.;
          rpa_force_dG_ij[i_coord][i_atom] =0.;
          rpa_force_dG_lk[i_coord][i_atom] =0.;
          rpa_force_dG_kl[i_coord][i_atom] =0.;
          rpa_force_dG_il[i_coord][i_atom] =0.;
          rpa_force_dG_li[i_coord][i_atom] =0.;
    }
    }   

        compute_self_energy_ji_times_dG(wc_times_v_R_tau, Cs, R_period, is, Mu, Nu, tau, R, itau, iR, rpa_force_dG_ji);
        compute_self_energy_jk_times_dG(wc_times_v_R_tau, Cs, R_period, is, Mu, Nu, tau, R, itau, iR, rpa_force_dG_kj);
        compute_self_energy_lk_times_dG(wc_times_v_R_tau, Cs, R_period, is, Mu, Nu, tau, R, itau, iR, rpa_force_dG_lk);
        compute_self_energy_li_times_dG(wc_times_v_R_tau, Cs, R_period, is, Mu, Nu, tau, R, itau, iR, rpa_force_dG_li);

	if (Mu != Nu)
        {		
        compute_self_energy_ji_times_dG(wc_times_v_R_tau, Cs, R_period, is, Nu, Mu, tau, R, itau, iR, rpa_force_dG_ij);
        compute_self_energy_jk_times_dG(wc_times_v_R_tau, Cs, R_period, is, Nu, Mu, tau, R, itau, iR, rpa_force_dG_jk);
        compute_self_energy_lk_times_dG(wc_times_v_R_tau, Cs, R_period, is, Nu, Mu, tau, R, itau, iR, rpa_force_dG_kl);
        compute_self_energy_li_times_dG(wc_times_v_R_tau, Cs, R_period, is, Nu, Mu, tau, R, itau, iR, rpa_force_dG_il);
	}



   for (int i_atom =0; i_atom !=natom; i_atom ++ )
   {		   
   for (int i_coord=1; i_coord<=3; i_coord++)
   {
    rpa_force_dG[i_coord][i_atom] = rpa_force_dG[i_coord][i_atom] +  rpa_force_dG_kj[i_coord][i_atom] + rpa_force_dG_jk[i_coord][i_atom] + 
                                    rpa_force_dG_ji[i_coord][i_atom] + rpa_force_dG_ij[i_coord][i_atom] + 
                                    rpa_force_dG_lk[i_coord][i_atom] + rpa_force_dG_kl[i_coord][i_atom] +
                                    rpa_force_dG_li[i_coord][i_atom] + rpa_force_dG_il[i_coord][i_atom]; 
    }
    }   


    }
      }
     }


         cout << "rpa force dG is " << endl;

      for (auto i_atom=0; i_atom<natom; i_atom++)
      {	
      cout << "atom   "  << i_atom << "  "  <<rpa_force_dG[1][i_atom]*constant_1 << "  " << rpa_force_dG[2][i_atom]*constant_1 << "  " << rpa_force_dG[3][i_atom]*constant_1 << endl; 
      for (int i_coord=1; i_coord<=3; i_coord++)
      {
       rpa_force[i_coord][i_atom+1] +=(-rpa_force_dG[i_coord][i_atom]);
      }
      }


	cout << "dV + dC + dG force" << endl;

        for (auto i_atom=1; i_atom<=natom; i_atom++)
	{	
        cout << "atom   "  << i_atom << "  "  <<rpa_force[1][i_atom]*constant_1 << "  " << rpa_force[2][i_atom]*constant_1 << "  " << rpa_force[3][i_atom]*constant_1 << endl; 
	}


   // ***********************************************************************************************************     

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


void  d_Chi0::find_n_homo_k_n_lumo_k(const MeanField &mf)
 {
    const auto nspins = mf.get_n_spins();
    double   myspin_double=static_cast<double>(nspins);	   
    const auto nkpts = mf.get_n_kpoints();
    const auto nbands = mf.get_n_bands();
    n_homo_k.create(nkpts, nspins);
    n_lumo_k.create(nkpts, nspins);
    n_homo_k.zero_out();
    n_lumo_k.zero_out();

    for (int is = 0; is != nspins; is++)
    {
        auto wg = mf.get_weight()[is];
            wg *=nkpts;
    for (int ik = 0; ik != nkpts; ik++)
    {
    for ( int ib = 0; ib != nbands; ib++)
    {
    if (wg(ik, ib) > 0.00000001)
    {
    n_homo_k(ik, is)=ib;
    }	    
    }

    for ( int ib =nbands-1; ib >= 0; ib--)
    {
    if (wg(ik, ib) <  2.0/myspin_double)
    {
     n_lumo_k(ik, is)=ib;
    }
    }
    
    }  
    }

 }


void d_Chi0::build_d_gf_full_unocc_R_plus_tau(Vector3_Order<int> R, double tau, const MeanField &mf)
{
    prof.start("cal_d_Green_func");
    const auto nkpts = mf.get_n_kpoints();
    const auto nspins = mf.get_n_spins();
    const auto nbands = mf.get_n_bands();
    const auto naos = mf.get_n_aos();
    assert (tau != 0);
    d_gf_save = d_gf_discard = 0;
    d_gf_save_x = d_gf_discard_x = 0;
    d_gf_save_y = d_gf_discard_y = 0;
    d_gf_save_z = d_gf_discard_z = 0;
    // temporary Green's function
    matrix gf_Rt_is_global(naos, naos);
    map<int, map<int, matrix>> d_gf_Rt_is_global; //(naos, naos);

    for (int is = 0; is != nspins; is++)
    {
        gf_Rt_is_global.zero_out();
        auto wg = mf.get_weight()[is];
            wg *=0.5 *nspins;
        for (int ik = 0; ik != nkpts; ik++)
            {
        for ( int ib = n_lumo_k(ik,is); ib != nbands; ib++)
            {
                wg(ik, ib) = 1.0 / nkpts  - wg(ik, ib);//
                if (wg(ik, ib) < 0) wg(ik, ib)= 0;
	    }
            }

        matrix scale(nkpts, nbands);
	scale.zero_out();
        map<int, map<int, matrix>>  d_scale;

        for (int i_coord=1; i_coord<=3; i_coord++) 
        {
        for (int i_atom=1; i_atom<=natom; i_atom++) 
	{
        d_scale[i_coord][i_atom].create(nkpts, nbands); 
        d_scale[i_coord][i_atom].zero_out(); 
	d_gf_Rt_is_global[i_coord][i_atom].create(naos, naos);
	d_gf_Rt_is_global[i_coord][i_atom].zero_out(); 
	}
        }		


        for (int ik = 0; ik != nkpts; ik++)
            {
        for ( int ib = n_lumo_k(ik,is); ib != nbands; ib++)
            {
           scale(ik, ib) = - 0.5 * tau * (mf.get_eigenvals()[is](ik, ib) - mf.get_efermi());
            // NOTE: enforce non-positive phase
            if ( scale(ik,ib) > 0) scale(ik,ib) = 0;
            scale(ik,ib) = std::exp(scale(ik,ib)) * wg(ik,ib);


        for (int i_coord=1; i_coord<=3; i_coord++) 
        {
        for (int i_atom=1; i_atom<=natom; i_atom++) 
	{
        d_scale[i_coord][i_atom](ik,ib)= - 0.5 * tau * scale(ik,ib) * meanfield.get_d_eigenvals()[i_coord][i_atom](ik,ib); 
        }
        }

	}
        }

        for (int ik = 0; ik != nkpts; ik++)
        {
            double ang = - klist[ik] * (R * latvec) * TWO_PI;
            complex<double> kphase = complex<double>(cos(ang), sin(ang));
            auto scaled_wfc_conj = conj(mf.get_eigenvectors()[is][ik]);
            for ( int ib = 0; ib != nbands; ib++)
                LapackConnector::scal(naos, scale(ik, ib), scaled_wfc_conj.c + naos * ib, 1);
            gf_Rt_is_global += (kphase * transpose(mf.get_eigenvectors()[is][ik], false) * scaled_wfc_conj).real();

        for (int i_coord=1; i_coord<=3; i_coord++) 
        {
        for (int i_atom=1; i_atom<=natom; i_atom++) 
	{
          auto scaled_wfc_conj1 = conj(mf.get_eigenvectors()[is][ik]);
          auto scaled_wfc_conj2 = conj(meanfield.get_d_eigenvectors()[i_coord-1][i_atom-1][ik]);
          for ( int ib = 0; ib != nbands; ib++)
	  {
          LapackConnector::scal(naos, d_scale[i_coord][i_atom](ik, ib), scaled_wfc_conj1.c + naos * ib, 1);
          LapackConnector::scal(naos, scale(ik, ib), scaled_wfc_conj2.c + naos * ib, 1);
	  }
          d_gf_Rt_is_global[i_coord][i_atom] += (kphase * transpose(mf.get_eigenvectors()[is][ik], false) * scaled_wfc_conj1).real();
          d_gf_Rt_is_global[i_coord][i_atom] += (kphase * transpose(mf.get_eigenvectors()[is][ik], false) * scaled_wfc_conj2).real();
          d_gf_Rt_is_global[i_coord][i_atom] += (kphase * transpose(meanfield.get_d_eigenvectors()[i_coord-1][i_atom-1][ik], false) *
			                         scaled_wfc_conj).real();
        }
	}

        }
//     cout << "bundesha is greate" << endl; 
       // if ( tau < 0 ) gf_Rt_is_global *= -1.;
        omp_lock_t gf_lock;
        omp_init_lock(&gf_lock);  
#pragma omp parallel for schedule(dynamic)
        for ( int I = 0; I != natom; I++)
        {
            const auto I_num = atom_nw[I];
            for (int J = 0; J != natom; J++)
            {
                const size_t J_num = atom_nw[J];
                matrix tmp_green(I_num, J_num);

                map<int, map<int, matrix>>  d_tmp_green;
            for ( int K = 0; K != natom; K++)
	    {
		d_tmp_green[1][K].create(I_num, J_num);
		d_tmp_green[2][K].create(I_num, J_num);
		d_tmp_green[3][K].create(I_num, J_num);
	    }

                for (size_t i = 0; i != I_num; i++)
                {
                    size_t i_glo = atom_iw_loc2glo(I, i);
                    for (size_t j = 0; j != J_num; j++)
                    {
                        size_t j_glo = atom_iw_loc2glo(J, j);
                        tmp_green(i, j) = gf_Rt_is_global(i_glo, j_glo);
               for ( int K = 0; K != natom; K++)
                {
                d_tmp_green[1][K](i, j) = d_gf_Rt_is_global[1][K+1](i_glo, j_glo);
                d_tmp_green[2][K](i, j) = d_gf_Rt_is_global[2][K+1](i_glo, j_glo);
                d_tmp_green[3][K](i, j) = d_gf_Rt_is_global[3][K+1](i_glo, j_glo);
                }
                    }
                }
                if (tmp_green.absmax() > gf_R_threshold)
                {
                    omp_set_lock(&gf_lock);
                    gf_is_R_tau[is][I][J][R][tau] = std::move(tmp_green);
                    omp_unset_lock(&gf_lock);
                    d_gf_save++;
                }
                else
                {
                    d_gf_discard++;
                }
             
               for ( int K = 0; K != natom; K++)
	       {
                if (d_tmp_green[1][K].absmax() > gf_R_threshold)
		{
                    d_gf_R_tau_full[1][K][I][J][R][tau] = std::move(d_tmp_green[1][K]);
                    d_gf_save_x++;
		}
                else
		{
		d_gf_discard_x++;	
		}		
                if (d_tmp_green[2][K].absmax() > gf_R_threshold)
		{
                    d_gf_R_tau_full[2][K][I][J][R][tau] = std::move(d_tmp_green[2][K]);
                    d_gf_save_y++;
		}
                else
		{
		d_gf_discard_y++;	
		}		

                if (d_tmp_green[3][K].absmax() > gf_R_threshold)
		{
                    d_gf_R_tau_full[3][K][I][J][R][tau] = std::move(d_tmp_green[3][K]);
                    d_gf_save_z++;
		}
                else
		{
		d_gf_discard_z++;	
		}		
                                 
               }


            }
        }
        omp_destroy_lock(&gf_lock);
#pragma omp barrier
    }
    prof.stop("cal_d_Green_func");
}


void d_Chi0::build_d_gf_full_occ_R_minus_tau(Vector3_Order<int> R, double tau, const MeanField &mf)
{
    prof.start("cal_d_Green_func");
    const auto nkpts = mf.get_n_kpoints();
    const auto nspins = mf.get_n_spins();
    const auto nbands = mf.get_n_bands();
    const auto naos = mf.get_n_aos();
    assert (tau != 0);
    d_gf_save = d_gf_discard = 0;
    d_gf_save_x = d_gf_discard_x = 0;
    d_gf_save_y = d_gf_discard_y = 0;
    d_gf_save_z = d_gf_discard_z = 0;
    // temporary Green's function
    matrix gf_Rt_is_global(naos, naos);
    map<int, map<int, matrix>> d_gf_Rt_is_global; //(naos, naos);

    for (int is = 0; is != nspins; is++)
    {
        gf_Rt_is_global.zero_out();
        auto wg = mf.get_weight()[is];
            wg *=0.5 *nspins;

        matrix scale(nkpts, nbands);
	scale.zero_out();
        map<int, map<int, matrix>>  d_scale;

        for (int i_coord=1; i_coord<=3; i_coord++) 
        {
        for (int i_atom=1; i_atom<=natom; i_atom++) 
	{
        d_scale[i_coord][i_atom].create(nkpts, nbands); 
        d_scale[i_coord][i_atom].zero_out(); 
	d_gf_Rt_is_global[i_coord][i_atom].create(naos, naos); 
	d_gf_Rt_is_global[i_coord][i_atom].zero_out(); 
	}
        }		


        for (int ik = 0; ik != nkpts; ik++)
            {
        for ( int ib = 0; ib <= n_homo_k(ik,is); ib++)
            {
           scale(ik, ib) = - 0.5 * tau * (mf.get_eigenvals()[is](ik, ib) - mf.get_efermi());
            // NOTE: enforce non-positive phase
            if ( scale(ik,ib) > 0) scale(ik,ib) = 0;
            scale(ik,ib) = std::exp(scale(ik,ib)) * wg(ik,ib);


        for (int i_coord=1; i_coord<=3; i_coord++) 
        {
        for (int i_atom=1; i_atom<=natom; i_atom++) 
	{
        d_scale[i_coord][i_atom](ik,ib)= - 0.5 * tau * scale(ik,ib) * meanfield.get_d_eigenvals()[i_coord][i_atom](ik,ib); 
        }
        }

	}
        }

        for (int ik = 0; ik != nkpts; ik++)
        {
            double ang = - klist[ik] * (R * latvec) * TWO_PI;
            complex<double> kphase = complex<double>(cos(ang), sin(ang));
            auto scaled_wfc_conj = conj(mf.get_eigenvectors()[is][ik]);
            for ( int ib=0; ib != nbands; ib++)
                LapackConnector::scal(naos, scale(ik, ib), scaled_wfc_conj.c + naos * ib, 1);
            gf_Rt_is_global += (kphase * transpose(mf.get_eigenvectors()[is][ik], false) * scaled_wfc_conj).real();

        for (int i_coord=1; i_coord<=3; i_coord++) 
        {
        for (int i_atom=1; i_atom<=natom; i_atom++) 
	{
         auto scaled_wfc_conj1 = conj(mf.get_eigenvectors()[is][ik]);
         auto scaled_wfc_conj2 = conj(meanfield.get_d_eigenvectors()[i_coord-1][i_atom-1][ik]);
         for ( int ib = 0; ib != nbands;  ib++)
	 {
         LapackConnector::scal(naos, d_scale[i_coord][i_atom](ik, ib), scaled_wfc_conj1.c + naos * ib, 1);
         LapackConnector::scal(naos, scale(ik, ib), scaled_wfc_conj2.c + naos * ib, 1);
	 }
         d_gf_Rt_is_global[i_coord][i_atom] -= (kphase * transpose(mf.get_eigenvectors()[is][ik], false) * scaled_wfc_conj1).real();
         d_gf_Rt_is_global[i_coord][i_atom] -= (kphase * transpose(mf.get_eigenvectors()[is][ik], false) * scaled_wfc_conj2).real();
         d_gf_Rt_is_global[i_coord][i_atom] -= (kphase * transpose(meanfield.get_d_eigenvectors()[i_coord-1][i_atom-1][ik], false) *
			                         scaled_wfc_conj).real();
        }
	}

        }
//     cout << "bundesha is greate" << endl; 
         gf_Rt_is_global *= -1.;
//         d_gf_Rt_is_global *= -1.;
        omp_lock_t gf_lock;
        omp_init_lock(&gf_lock);  
#pragma omp parallel for schedule(dynamic)
        for ( int I = 0; I != natom; I++)
        {
            const auto I_num = atom_nw[I];
            for (int J = 0; J != natom; J++)
            {
                const size_t J_num = atom_nw[J];
                matrix tmp_green(I_num, J_num);

                map<int, map<int, matrix>>  d_tmp_green;
            for (int K = 0; K != natom; K++)
	    {
		d_tmp_green[1][K].create(I_num, J_num);
		d_tmp_green[2][K].create(I_num, J_num);
		d_tmp_green[3][K].create(I_num, J_num);
	    }

                for (size_t i = 0; i != I_num; i++)
                {
                    size_t i_glo = atom_iw_loc2glo(I, i);
                    for (size_t j = 0; j != J_num; j++)
                    {
                        size_t j_glo = atom_iw_loc2glo(J, j);
                        tmp_green(i, j) = gf_Rt_is_global(i_glo, j_glo);
                  for (int K = 0; K != natom; K++)
		  {
                        d_tmp_green[1][K](i, j) = d_gf_Rt_is_global[1][K+1](i_glo, j_glo);
                        d_tmp_green[2][K](i, j) = d_gf_Rt_is_global[2][K+1](i_glo, j_glo);
                        d_tmp_green[3][K](i, j) = d_gf_Rt_is_global[3][K+1](i_glo, j_glo);
		  }
                    }
                }
                if (tmp_green.absmax() > gf_R_threshold)
                {
                    omp_set_lock(&gf_lock);
                    gf_is_R_tau[is][I][J][R][tau] = std::move(tmp_green);
                    omp_unset_lock(&gf_lock);
                    d_gf_save++;
                }
                else
                {
                    d_gf_discard++;
                }

                  for (int K = 0; K != natom; K++)
		  {
                if (d_tmp_green[1][K].absmax() > gf_R_threshold)
		{
                    d_gf_R_tau_full[1][K][I][J][R][tau] = std::move(d_tmp_green[1][K]);
                    d_gf_save_x++;
		}
		else
		{
	          d_gf_discard_x++;		
		}

                if (d_tmp_green[2][K].absmax() > gf_R_threshold)
		{
                    d_gf_R_tau_full[2][K][I][J][R][tau] = std::move(d_tmp_green[2][K]);
                    d_gf_save_y++;
		}
		else
		{
	          d_gf_discard_y++;		
		}

                if (d_tmp_green[3][K].absmax() > gf_R_threshold)
		{
                    d_gf_R_tau_full[3][K][I][J][R][tau] = std::move(d_tmp_green[3][K]);
                    d_gf_save_z++;
		}
		else
		{
	          d_gf_discard_z++;		
		}

		  }


            }
        }
        omp_destroy_lock(&gf_lock);
#pragma omp barrier
    }
    prof.stop("cal_d_Green_func");
}


void  d_Chi0::compute_self_energy_jk_times_dG(atom_mapping< map<int, map<int, ComplexMatrix>>>::pair_t_old & wc_times_v_R_tau,
	                                    const atpair_R_mat_t &LRI_Cs,
                                            const Vector3_Order<int> &R_period, 
                                            int spin_channel, atom_t Mu, atom_t Nu, double tau, 
				            Vector3_Order<int> R, int itau, int iR, 
					    map<int, map<int, double>>  & rpa_force_dG)
   {

    assert ( tau > 0 );
    // Local RI requires
    const atom_t I_index = Mu;
    const atom_t J_index = Nu;

    const size_t i_num = atom_nw[Mu];
    const size_t mu_num = atom_mu[Mu];

    const size_t j_num = atom_nw[Nu];
    const size_t nu_num = atom_mu[Nu];


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

       for (int i_coord=1; i_coord<=3; i_coord++)
       {
       for (int i_atom =0; i_atom !=natom; i_atom ++ )
       {		   
       rpa_force_dG[i_coord][i_atom] = 0.;
       }
       }

   

    for (const auto &K_pair : Cs[I_index])
    {
        const auto K_index = K_pair.first;
        const size_t k_num = atom_nw[K_index];

	ComplexMatrix    Wc_times_Cs(nu_num, i_num*k_num);  
       for (const auto &R1_index : K_pair.second)
        {
           const auto R1 = R1_index.first;
           const auto &Cs_mat1 = R1_index.second;
          Vector3_Order<int> R_temp_2(Vector3_Order<int>(R - R1) % R_period);
        // Wc_times_Cs(nu, i*k) 
        Wc_times_Cs  = wc_times_v_R_tau[J_index][I_index][iR][itau] * ComplexMatrix(transpose(*Cs_mat1)) ; 
	// Wc_times_Cs_trans(i_num*k_num, nu_num) 
	auto  Wc_times_Cs_trans=transpose(Wc_times_Cs, false);
//        ComplexMatrix reshap_mat11(reshape_complxmat_12(i_num, k_num, nu_num, Wc_times_Cs_trans)); // (n1*n2,n3) ---> (n1,n2*n3)
//        ComplexMatrix Wc_times_Cs_rs(reshape_complxmat(i_num, k_num, nu_num, reshap_mat11)); //(n1,n2*n3) -> (n2,n1*n3)
        //  (i*k, nu) ---> (k, i*nu) 
        ComplexMatrix   Wc_times_Cs_rs(reshape_complxmat_13(i_num, k_num, nu_num, Wc_times_Cs_trans));




      matrix Z(i_num, j_num*nu_num); 
      matrix Z_conj(i_num, j_num*nu_num); 
	Z.zero_out();
	Z_conj.zero_out();

        ComplexMatrix  self_energy_jk_plus_tau(j_num, k_num); 
        ComplexMatrix  self_energy_jk_minus_tau(j_num, k_num); 


      for (const auto & L_pair : Cs[J_index] ) 
      {
        const auto L_index =L_pair.first;
        const size_t l_num = atom_nw[L_index];
      if (gf_R_tau.at(I_index).count(L_index))
       {		
       for (const auto &R2_index : L_pair.second)
       {
           const auto R2 = R2_index.first;
           const auto &Cs_mat2 = R2_index.second;
           Vector3_Order<int> R_temp_3(Vector3_Order<int>(R + R2) % R_period);
       if ( gf_R_tau.at(I_index).at(L_index).count(R_temp_3))
       {
         assert(j_num * l_num == (*Cs_mat2).nr);
        matrix Cs2_reshape(reshape_Cs(j_num, l_num, nu_num, Cs_mat2));


           if (gf_R_tau.at(I_index).at(L_index).at(R_temp_3).count(tau))
           {
           // Z(i, j*nu) = G(i, l) * Cs(l, j*nu) 
               Z += gf_R_tau.at(I_index).at(L_index).at(R_temp_3).at(tau)* Cs2_reshape;
           }

           if (gf_R_tau.at(I_index).at(L_index).at(R_temp_3).count(-tau))
           {
	   // Z_conj(i, j*nu) = G(i, l) * Cs(l, j*nu) 
               Z_conj += gf_R_tau.at(I_index).at(L_index).at(R_temp_3).at(-tau) * Cs2_reshape;
           }

       }

       }
      }
      }

          // Z(i, j*nu) ---> (j, i*nu) 
          matrix Z_rs(reshape_mat(i_num, j_num, nu_num, Z));                
          matrix Z_conj_rs(reshape_mat(i_num, j_num, nu_num, Z_conj));
           
          self_energy_jk_plus_tau.zero_out();
          self_energy_jk_minus_tau.zero_out();
              //(j, k) = Z(j, i*nu) * Wc_times_Cs_rs(k, i*nu)   
         self_energy_jk_plus_tau = ComplexMatrix(Z_rs) * transpose(Wc_times_Cs_rs, false); 
         self_energy_jk_minus_tau = ComplexMatrix(Z_conj_rs) * transpose(Wc_times_Cs_rs, false); 

         ComplexMatrix self_ene_times_dG(j_num, j_num);
       for (int i_coord=1; i_coord<=3; i_coord++)
       {
       for (int i_atom =0; i_atom !=natom; i_atom ++ )
       {		   
	 self_ene_times_dG.zero_out();
       if (d_gf_R_tau_full.at(i_coord).at(i_atom).at(K_index).at(J_index).at(R_temp_2).count(tau))
       {
        self_ene_times_dG = self_energy_jk_minus_tau * ComplexMatrix(d_gf_R_tau_full[i_coord][i_atom][K_index][J_index][R_temp_2][tau]);
       }
       if (d_gf_R_tau_full.at(i_coord).at(i_atom).at(K_index).at(J_index).at(R_temp_2).count(-tau))
       {
        self_ene_times_dG += self_energy_jk_plus_tau * ComplexMatrix(d_gf_R_tau_full[i_coord][i_atom][K_index][J_index][R_temp_2][-tau]);
       }

        complex<double> trace_self_mat_dG=trace(self_ene_times_dG); 
        rpa_force_dG[i_coord][i_atom] +=  2.*trace_self_mat_dG.real();	    
       }
       }
       

       }

     }


   }


void  d_Chi0::compute_self_energy_lk_times_dG(atom_mapping< map<int, map<int, ComplexMatrix>>>::pair_t_old & wc_times_v_R_tau,
	                                    const atpair_R_mat_t &LRI_Cs,
                                            const Vector3_Order<int> &R_period, 
                                            int spin_channel, atom_t Mu, atom_t Nu, double tau, 
				            Vector3_Order<int> R, int itau, int iR, 
					    map<int, map<int, double>>  & rpa_force_dG)
   {

    assert ( tau > 0 );
    // Local RI requires
    const atom_t I_index = Mu;
    const atom_t J_index = Nu;

    const size_t i_num = atom_nw[Mu];
    const size_t mu_num = atom_mu[Mu];

    const size_t j_num = atom_nw[Nu];
    const size_t nu_num = atom_mu[Nu];


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

       for (int i_coord=1; i_coord<=3; i_coord++)
       {
       for (int i_atom =0; i_atom !=natom; i_atom ++ )
       {		   
       rpa_force_dG[i_coord][i_atom] = 0.;
       }
       }
   

    for (const auto &K_pair : Cs[I_index])
    {
        const auto K_index = K_pair.first;
        const size_t k_num = atom_nw[K_index];

	ComplexMatrix    Wc_times_Cs(nu_num, i_num*k_num);  
       for (const auto &R1_index : K_pair.second)
        {
           const auto R1 = R1_index.first;
           const auto &Cs_mat1 = R1_index.second;
          Vector3_Order<int> R_temp_2(Vector3_Order<int>(R - R1) % R_period);
        // Wc_times_Cs(nu, i*k) 
        Wc_times_Cs  = wc_times_v_R_tau[J_index][I_index][iR][itau] * ComplexMatrix(transpose(*Cs_mat1)) ; 
	// Wc_times_Cs_trans(i_num*k_num, nu_num) 
	auto  Wc_times_Cs_trans=transpose(Wc_times_Cs, false);
        //  (i*k, nu) ---> (k, i*nu) 
        ComplexMatrix   Wc_times_Cs_rs(reshape_complxmat_13(i_num, k_num, nu_num, Wc_times_Cs_trans));


       if (flag_G_IJRt || flag_G_IJRNt ) 
       {
      for (const auto & L_pair : Cs[J_index] ) 
      {
        const auto L_index =L_pair.first;
        const size_t l_num = atom_nw[L_index];
        ComplexMatrix  self_energy_lk_plus_tau(l_num, k_num); 
        ComplexMatrix  self_energy_lk_minus_tau(l_num, k_num); 
        matrix M(i_num, l_num*nu_num); 
        matrix M_conj(i_num, l_num*nu_num); 
	M.zero_out();
	M_conj.zero_out();
       for (const auto &R2_index : L_pair.second)
       {
           const auto R2 = R2_index.first;
           const auto &Cs_mat2 = R2_index.second;
           Vector3_Order<int> R_temp_3(Vector3_Order<int>(R + R2) % R_period);
         assert(j_num * l_num == (*Cs_mat2).nr);
       // (j*l, nu) ---> (j, l*mu) 
        matrix Cs2_reshape(reshape_Cs_12(j_num, l_num, nu_num, Cs_mat2));

           if (flag_G_IJRt)
           {
           // M(i, l*nu) = G(i, j) * Cs(j, l*nu) 
               M = gf_R_tau.at(I_index).at(J_index).at(R).at(tau)* Cs2_reshape;
           }
           if (flag_G_IJRNt)
           {
	   // M_conj(i, l*nu) = G(i, j) * Cs(j, l*nu) 
               M_conj = gf_R_tau.at(I_index).at(J_index).at(R).at(-tau) * Cs2_reshape;
           }

	   // M_rs(l, i*nu)  
         matrix    M_rs(reshape_mat(i_num, l_num, nu_num, M));
         matrix    M_conj_rs(reshape_mat(i_num, l_num, nu_num, M_conj));
          
         self_energy_lk_plus_tau.zero_out();
	 self_energy_lk_minus_tau.zero_out();

	 self_energy_lk_plus_tau = ComplexMatrix(M_rs) * transpose(Wc_times_Cs_rs, false); 
	 self_energy_lk_minus_tau = ComplexMatrix(M_conj_rs) * transpose(Wc_times_Cs_rs, false); 

         Vector3_Order<int> R_temp_1(Vector3_Order<int>(R + R2 - R1) % R_period);

         ComplexMatrix self_ene_times_dG(l_num, l_num);

       for (int i_coord=1; i_coord<=3; i_coord++)
       {
       for (int i_atom =0; i_atom !=natom; i_atom ++ )
       {		   
	   self_ene_times_dG.zero_out();
       if (d_gf_R_tau_full.at(i_coord).at(i_atom).at(K_index).at(L_index).at(R_temp_1).count(tau)) 
        {
         self_ene_times_dG = self_energy_lk_minus_tau *  ComplexMatrix(d_gf_R_tau_full[i_coord][i_atom][K_index][L_index][R_temp_1][tau]);
	}
       if (d_gf_R_tau_full.at(i_coord).at(i_atom).at(K_index).at(L_index).at(R_temp_1).count(-tau)) 
       {
         self_ene_times_dG += self_energy_lk_plus_tau *  ComplexMatrix(d_gf_R_tau_full[i_coord][i_atom][K_index][L_index][R_temp_1][-tau]);
       }
        complex<double> trace_self_mat_dG=trace(self_ene_times_dG); 
        rpa_force_dG[i_coord][i_atom] +=  2.*trace_self_mat_dG.real();	    
       }
       }


       }
       }
       }
       

       }
     }


   }



void  d_Chi0::compute_self_energy_li_times_dG(atom_mapping< map<int, map<int, ComplexMatrix>>>::pair_t_old & wc_times_v_R_tau,
	                                    const atpair_R_mat_t &LRI_Cs,
                                            const Vector3_Order<int> &R_period, 
                                            int spin_channel, atom_t Mu, atom_t Nu, double tau, 
				            Vector3_Order<int> R, int itau, int iR, 
                                            map<int, map<int, double>> & rpa_force_dG)


   {

    assert ( tau > 0 );
    // Local RI requires
    const atom_t I_index = Mu;
    const atom_t J_index = Nu;

    const size_t i_num = atom_nw[Mu];
    const size_t mu_num = atom_mu[Mu];

    const size_t j_num = atom_nw[Nu];
    const size_t nu_num = atom_mu[Nu];


    int flag_G_IJRt = 0;
    int flag_G_IJRNt = 0;

    const auto & gf_R_tau = gf_is_R_tau.at(spin_channel);
//    if (gf_R_tau.at(I_index).count(J_index))
//        if (gf_R_tau.at(I_index).at(J_index).count(R))
//        {
//            if (gf_R_tau.at(I_index).at(J_index).at(R).count(tau))
//                flag_G_IJRt = 1;
//            if (gf_R_tau.at(I_index).at(J_index).at(R).count(-tau))
//                flag_G_IJRNt = 1;
//        }



    for (const auto &K_pair : Cs[I_index])
    {
        const auto K_index = K_pair.first;
        const size_t k_num = atom_nw[K_index];
	ComplexMatrix    Wc_times_Cs(nu_num, i_num*k_num);  
     if (gf_R_tau.at(K_index).count(J_index))
      {	     
       for (const auto &R1_index : K_pair.second)
        {
           const auto R1 = R1_index.first;
           const auto &Cs_mat1 = R1_index.second;
           Vector3_Order<int> R_temp_2(Vector3_Order<int>(R - R1) % R_period);
       if (gf_R_tau.at(K_index).at(J_index).count(R_temp_2))
       {	       
        // Wc_times_Cs(nu, i*k) 
        Wc_times_Cs  = wc_times_v_R_tau[J_index][I_index][iR][itau] * ComplexMatrix(transpose(*Cs_mat1)) ; 
	// Wc_times_Cs_trans(i_num*k_num, nu_num) 
	auto  Wc_times_Cs_trans=transpose(Wc_times_Cs, false);
        //  (i*k, nu) ---> (i, k*nu) 
        ComplexMatrix   Wc_times_Cs_rs(reshape_complxmat_12(i_num, k_num, nu_num, Wc_times_Cs_trans));


      for (const auto & L_pair : Cs[J_index] ) 
      {
        const auto L_index =L_pair.first;
        const size_t l_num = atom_nw[L_index];

        ComplexMatrix self_energy_li_plus_tau(l_num, i_num);
        ComplexMatrix self_energy_li_minus_tau(l_num, i_num);

        matrix X_R2(k_num, l_num * nu_num);
        matrix X_conj_R2(k_num, l_num * nu_num);

       for (const auto &R2_index : L_pair.second)
       {
           const auto R2 = R2_index.first;
           const auto &Cs_mat2 = R2_index.second;
           Vector3_Order<int> R_temp_3(Vector3_Order<int>(R + R2) % R_period);
         assert(j_num * l_num == (*Cs_mat2).nr);
	 //(j*l, nu) ---> (j, l*nu) 
        matrix Cs2_reshape(reshape_Cs_12(j_num, l_num, nu_num, Cs_mat2));

       if (gf_R_tau.at(K_index).at(J_index).at(R_temp_2).count(tau))
       {   
           //X(k, l*\nu) 
           X_R2 = gf_R_tau.at(K_index).at(J_index).at(R_temp_2).at(tau) * Cs2_reshape;
       }
       if (gf_R_tau.at(K_index).at(J_index).at(R_temp_2).count(-tau))
       {
           //X(k, l*\nu) 
           X_conj_R2 = gf_R_tau.at(K_index).at(J_index).at(R_temp_2).at(-tau) * Cs2_reshape;
       }

	 //N_rs(l, k*\nu)
          matrix X_conj_R2_rs(reshape_mat(k_num, l_num, nu_num, X_conj_R2));
	 //N_rs(l, k*\nu)
          matrix X_R2_rs(reshape_mat(k_num, l_num, nu_num, X_R2));

          self_energy_li_plus_tau.zero_out();
          self_energy_li_minus_tau.zero_out();

	  self_energy_li_plus_tau  =  ComplexMatrix(X_R2_rs) * transpose(Wc_times_Cs_rs, false); 
	  self_energy_li_minus_tau =  ComplexMatrix(X_conj_R2_rs) * transpose(Wc_times_Cs_rs, false); 
          


          ComplexMatrix self_energy_times_dG(l_num, l_num); 


         for (int i_atom =0; i_atom !=natom; i_atom ++ )
         {		   
         for (int i_coord=1; i_coord<=3; i_coord++)
         {
            self_energy_times_dG.zero_out();
         if (d_gf_R_tau_full.at(i_coord).at(i_atom).at(I_index).at(L_index).at(R_temp_3).count(tau))
          {
          self_energy_times_dG = self_energy_li_minus_tau * ComplexMatrix(d_gf_R_tau_full[i_coord][i_atom][I_index][L_index][R_temp_3][tau]); 
          }
         if (d_gf_R_tau_full.at(i_coord).at(i_atom).at(I_index).at(L_index).at(R_temp_3).count(-tau))
          {
          self_energy_times_dG += self_energy_li_plus_tau * ComplexMatrix(d_gf_R_tau_full[i_coord][i_atom][I_index][L_index][R_temp_3][-tau]);  
          }
                  complex<double> trace_self_mat_dG=trace(self_energy_times_dG); 
                  rpa_force_dG[i_coord][i_atom] +=   2.*trace_self_mat_dG.real();	    
          }
         } 


       }
      }


       }
       }
    }
    }


   }

void  d_Chi0::compute_self_energy_ji_times_dG(atom_mapping< map<int, map<int, ComplexMatrix>>>::pair_t_old & wc_times_v_R_tau,
	                                    const atpair_R_mat_t &LRI_Cs,
                                            const Vector3_Order<int> &R_period, 
                                            int spin_channel, atom_t Mu, atom_t Nu, double tau, 
				            Vector3_Order<int> R, int itau, int iR, 
                                            map<int, map<int, double>> & rpa_force_dG)


   {

    assert ( tau > 0 );
    // Local RI requires
    const atom_t I_index = Mu;
    const atom_t J_index = Nu;

    const size_t i_num = atom_nw[Mu];
    const size_t mu_num = atom_mu[Mu];

    const size_t j_num = atom_nw[Nu];
    const size_t nu_num = atom_mu[Nu];


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



      ComplexMatrix self_energy_ji_plus_tau(j_num, i_num);
      ComplexMatrix self_energy_ji_minus_tau(j_num, i_num);
      self_energy_ji_plus_tau.zero_out();
      self_energy_ji_minus_tau.zero_out();


    for (const auto &K_pair : Cs[I_index])
    {
        const auto K_index = K_pair.first;
        const size_t k_num = atom_nw[K_index];

	ComplexMatrix    Wc_times_Cs(nu_num, i_num*k_num);  
       for (const auto &R1_index : K_pair.second)
        {
           const auto R1 = R1_index.first;
           const auto &Cs_mat1 = R1_index.second;
        // Wc_times_Cs(nu, i*k) 
        Wc_times_Cs  = wc_times_v_R_tau[J_index][I_index][iR][itau] * ComplexMatrix(transpose(*Cs_mat1)) ; 
	// Wc_times_Cs_trans(i_num*k_num, nu_num) 
	auto  Wc_times_Cs_trans=transpose(Wc_times_Cs, false);
        //  (i*k, nu) ---> (i, k*nu) 
        ComplexMatrix   Wc_times_Cs_rs(reshape_complxmat_12(i_num, k_num, nu_num, Wc_times_Cs_trans));

        matrix N_R2(k_num, j_num * nu_num);
        matrix N_conj_R2(k_num, j_num * nu_num);

      for (const auto & L_pair : Cs[J_index] ) 
      {
        const auto L_index =L_pair.first;
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
               //N(k, j*\nu) 
           N_R2 += gf_R_tau.at(K_index).at(L_index).at(R_temp_1).at(tau) * Cs2_reshape;
       }
       if (flag_G_IJRt && gf_R_tau.at(K_index).at(L_index).at(R_temp_1).count(-tau))
       {
               //N(k, j*\nu) 
           N_conj_R2 += gf_R_tau.at(K_index).at(L_index).at(R_temp_1).at(-tau) * Cs2_reshape;
       }


       }


       }
      }
      }
	   //N_rs(j, k*\nu)
            matrix N_conj_R2_rs(reshape_mat(k_num, j_num, nu_num, N_conj_R2));
	   //N_rs(j, k*\nu)
            matrix N_R2_rs(reshape_mat(k_num, j_num, nu_num, N_R2));

           if (flag_G_IJRt)
           {  
	    // self_energy_ji_minus_tau (j,i)= N_conj_R2_rs(j, k*nu) * transpose(Wc_times_Cs_rs(i, k*nu))
	    self_energy_ji_minus_tau += ComplexMatrix(N_conj_R2_rs) * transpose(Wc_times_Cs_rs, false);  
           }
           if (flag_G_IJRNt)
           {
	    // self_energy_ji_plus_tau (j,i)= N_R2_rs(j, k*nu) * transpose(Wc_times_Cs_rs(i, k*nu))
	    self_energy_ji_plus_tau += ComplexMatrix(N_R2_rs) * transpose(Wc_times_Cs_rs, false);  
           }

       }

     }


          self_energy_ji_minus_tau = 2.*self_energy_ji_minus_tau; 
	  self_energy_ji_plus_tau  = 2.*self_energy_ji_plus_tau;




    ComplexMatrix self_energy_times_dG(j_num, j_num); 



   for (int i_atom =0; i_atom !=natom; i_atom ++ )
   {		   
   for (int i_coord=1; i_coord<=3; i_coord++)
   {
      self_energy_times_dG.zero_out();
   if (d_gf_R_tau_full.at(i_coord).at(i_atom).at(Mu).at(Nu).at(R).count(tau))
    {
    self_energy_times_dG = self_energy_ji_minus_tau * ComplexMatrix(d_gf_R_tau_full[i_coord][i_atom][Mu][Nu][R][tau]); 
    }
   if (d_gf_R_tau_full.at(i_coord).at(i_atom).at(Mu).at(Nu).at(R).count(-tau))
    {
    self_energy_times_dG += self_energy_ji_plus_tau * ComplexMatrix(d_gf_R_tau_full[i_coord][i_atom][Mu][Nu][R][-tau]);  
    }
	    complex<double> trace_self_mat_dG=trace(self_energy_times_dG); 
            rpa_force_dG[i_coord][i_atom] =   trace_self_mat_dG.real();	    
    }
   } 



   }
