
#include  <map>  
#include  <iostream>
#include "vector3_order.h"
#include "ri.h"
#include "timefreq.h"
#include "parallel_mpi.h"
#include <vector>
#include "meanfield.h"
//#include "chi0.h"


  extern  vector<Vector3_Order<int>> my_Rlist;
  extern  vector<Vector3_Order<double>> qlist;
  extern  vector<atpair_t> atpairs_ABF;
  extern  Vector3_Order<int> R_period;
  extern  map<int, atom_mapping<map<Vector3_Order<int>, map<double, matrix>>>::pair_t_old> gf_is_R_tau;
  extern  map<int, atom_mapping<map<Vector3_Order<int>, map<double, matrix>>>::pair_t_old> d_gf_R_tau;
  extern  map<int, map<int, atom_mapping<map<Vector3_Order<int>, map<double, matrix>>>::pair_t_old>> d_gf_R_tau_full;
  extern  map<int, map<int, double>> rpa_force;
  extern  matrix  n_homo_k, n_lumo_k;
 class d_Chi0
 {
    public:
    double gf_R_threshold;
   //void   integrate_freq_wc_times_v(); 
   void  integrate_freq_wc_times_v(map<double,map<Vector3_Order<double>, ComplexMatrix>> &   wc_times_v,
		   const TFGrids & tfg, map<int, map<Vector3_Order<double>, ComplexMatrix>> & wc_times_v_out); 

   //void  integrate_time_wc_times_v_R_tau(const Chi0 &chi0,map<int, map<Vector3_Order<double>, ComplexMatrix>> & wc_times_v_out, 
   void  integrate_time_wc_times_v_R_tau(map<int, map<Vector3_Order<double>, ComplexMatrix>> & wc_times_v_out, 
		   const TFGrids & tfg, const MeanField & mf);//, map<int, map<int, ComplexMatrix>> & wc_times_v_R_tau); 

   //void compute_wc_times_O();
   //void compute_wc_times_O(const Chi0 &chi0, map<int, map<int, ComplexMatrix>> & wc_times_v_R_tau, 
   void compute_wc_times_O(atom_mapping< map<int, map<int, ComplexMatrix>>>::pair_t_old  & wc_times_v_R_tau, 
		             vector<pair<int,int>> & itauiRs_local, const TFGrids & tfg, const MeanField &mf );


void build_d_gf_full_unocc_R_plus_tau(Vector3_Order<int> R, double tau, const MeanField &mf);
void build_d_gf_full_occ_R_minus_tau(Vector3_Order<int> R, double tau, const MeanField &mf);


void  find_n_homo_k_n_lumo_k(const MeanField &mf);


void  compute_self_energy_ji_times_dG(atom_mapping< map<int, map<int, ComplexMatrix>>>::pair_t_old & wc_times_v_R_tau,
	                                    const atpair_R_mat_t &LRI_Cs,
                                            const Vector3_Order<int> &R_period, 
                                            int spin_channel, atom_t Mu, atom_t Nu, double tau, 
				            Vector3_Order<int> R, int itau, int iR,
                                            map<int, map<int, double>> & rpa_force_dG);


void  compute_self_energy_lk_times_dG(atom_mapping< map<int, map<int, ComplexMatrix>>>::pair_t_old & wc_times_v_R_tau,
	                                    const atpair_R_mat_t &LRI_Cs,
                                            const Vector3_Order<int> &R_period, 
                                            int spin_channel, atom_t Mu, atom_t Nu, double tau, 
				            Vector3_Order<int> R, int itau, int iR, 
					    map<int, map<int, double>>  & rpa_force_dG);


void  compute_self_energy_li_times_dG(atom_mapping< map<int, map<int, ComplexMatrix>>>::pair_t_old & wc_times_v_R_tau,
	                                    const atpair_R_mat_t &LRI_Cs,
                                            const Vector3_Order<int> &R_period, 
                                            int spin_channel, atom_t Mu, atom_t Nu, double tau, 
				            Vector3_Order<int> R, int itau, int iR, 
                                            map<int, map<int, double>> & rpa_force_dG);


void  compute_self_energy_jk_times_dG(atom_mapping< map<int, map<int, ComplexMatrix>>>::pair_t_old & wc_times_v_R_tau,
	                                    const atpair_R_mat_t &LRI_Cs,
                                            const Vector3_Order<int> &R_period, 
                                            int spin_channel, atom_t Mu, atom_t Nu, double tau, 
				            Vector3_Order<int> R, int itau, int iR, 
					    map<int, map<int, double>>  & rpa_force_dG);


//matrix compute_O_s_nu_tau_R_R1( const Chi0 & chi0, const atpair_R_mat_t &LRI_Cs,
matrix compute_O_s_nu_tau_R_R1(const atpair_R_mat_t &LRI_Cs,
                                       const Vector3_Order<int> &R_period, 
                                       int spin_channel,
                                       atom_t Mu, atom_t Nu, atom_t Ku, double tau, 
				       Vector3_Order<int> R, Vector3_Order<int> R1);


matrix compute_O_bar_s_nu_tau_R_R2(const atpair_R_mat_t &LRI_Cs,
                                       const Vector3_Order<int> &R_period, 
                                       int spin_channel,
                                       atom_t Mu, atom_t Nu, atom_t Lu, double tau, 
				       Vector3_Order<int> R, Vector3_Order<int> R2);



    private:
     vector<Vector3_Order<int>> Rlist_gf;
     size_t d_gf_discard;
     size_t d_gf_discard_x;
     size_t d_gf_discard_y;
     size_t d_gf_discard_z;
     size_t d_gf_save;
     size_t d_gf_save_x;
     size_t d_gf_save_y;
     size_t d_gf_save_z;
    // MeanField  mf;
//     map<int, atom_mapping<map<Vector3_Order<int>, map<double, matrix>>>::pair_t_old> gf_is_R_tau;
  //   void build_gf_Rt(Vector3_Order<int> R, double tau);
 };


