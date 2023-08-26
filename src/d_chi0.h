
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
  extern  map<int, map<int, double>> rpa_force;
 class d_Chi0
 {
    public:
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
//     map<int, atom_mapping<map<Vector3_Order<int>, map<double, matrix>>>::pair_t_old> gf_is_R_tau;
  //   void build_gf_Rt(Vector3_Order<int> R, double tau);
 };


