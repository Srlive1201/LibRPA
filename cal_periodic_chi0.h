#ifndef CAL_PERIODIC_CHI0_H
#define CAL_PERIODIC_CHI0_H

#include "matrix.h"
#include "complexmatrix.h"
#include <vector>
#include <map>
#include <set>
#include <memory>
#include "vector3.h"
#include "vector3_order.h"
#include "input.h"
#include <cstring>
class Cal_Periodic_Chi0
{
    //double e_const=2.718281828;
    
    std::vector<ComplexMatrix> green_k;		// green_k[ik](iw1,iw2);
    std::map<size_t,map<double,map<Vector3_Order<int>,ComplexMatrix>>> Green_function;  //Green_fuction[is][time_tau][R](iw1,iw2);
    std::map<double,map<Vector3_Order<int>,map<size_t,map<size_t,ComplexMatrix>>>> chi0;
    std::map<double,map<Vector3_Order<int>,map<size_t,map<size_t,ComplexMatrix>>>> chi0_freq;
    std::map<double,map<Vector3_Order<double>,map<size_t,map<size_t,ComplexMatrix>>>> chi0_k;
    std::map<double,map<Vector3_Order<double>,map<size_t,map<size_t,ComplexMatrix>>>> pi_k; 
    //std::map<double,map<Abfs::Vector3_Order<double>,map<size_t,map<size_t,ComplexMatrix>>>> pi_k;
    //std::map<double,map<Abfs::Vector3_Order<int>,ComplexMatrix>> Green_out;
    int kpoints_num;
    std::set<size_t> iwt2iat_row;
    std::set<size_t> iwt2iat_col;
    
    size_t all_mu,all_nu;

    //map<double,double> freq_grid;
    //vector<Vector3<double>> R_grid;
    //vector<double> time_tau_grid;
    //map<size_t,map<size_t,map<Vector3_Order<int>,std::shared_ptr<matrix>>>> Cs_m;
    
    Vector3_Order<int> R_periodic;
    void init();
    void cal_Green_func_element_k( const matrix &wg, const double time_tau);
	void cal_Green_func_R_tau(const double &time_tau, const Vector3_Order<int> &R, const Vector3<double>* kvec_c);
    set<Vector3_Order<int>> construct_R_grid();
    map<double,double> construct_gauss_grid(const int Npoints);
    double get_band_gap();
    map<double,double> read_file_grid(const string &file_path, const char type);
    void cal_chi0_element(const double &time_tau, const Vector3_Order<int> &R); 
    void cal_chi0_element_origin(const double &time_tau, const Vector3_Order<int> &R);
    //void cal_local_Cs();
    matrix reshape_Cs(const size_t n1, const size_t n2, const size_t n3, const std::shared_ptr<matrix> &Cs);
    matrix reshape_dim_Cs(const size_t n1, const size_t n2, const size_t n3, const std::shared_ptr<matrix> &Cs);//(n1*n2,n3) -> (n1,n2*n3)
    ComplexMatrix reshape_complexmat(const size_t n1, const size_t n2, const size_t n3, const ComplexMatrix &mat);
    ComplexMatrix reshape_complexmat_21(const size_t n1, const size_t n2, const size_t n3, const ComplexMatrix &mat); //(n1,n2*n3) -> (n1*n2,n3)
    void generate_Green_part(const size_t &is, const double &time_tau, const Vector3_Order<int> &R, const size_t K_index, const size_t L_index, ComplexMatrix &Green_kl, ComplexMatrix &Green_kl_conj);
    Vector3_Order<int> periodic_R(const Vector3_Order<int> &R); 
//    void FT_to_chi0_freq(const set<Vector3_Order<int>> &R_grid, map<double,double> &time_grid, map<double,double> &freq_grid );
    void Cosine_to_chi0_freq(const set<Vector3_Order<int>> &R_grid, map<double,double> &time_grid, map<double,double> &freq_grid );
    map<double,map<Vector3_Order<double>,map<size_t,map<size_t,ComplexMatrix>>>> cal_pi_k_use_aims_vq();
    void FT_R_to_K_chi0();
    //void cal_pi_freq();
/*    
    std::map<double,map<Vector3_Order<double>,map<size_t,map<size_t,ComplexMatrix>>>> cal_pi_k();
    std::map<double,map<Vector3_Order<double>,map<size_t,map<size_t,ComplexMatrix>>>> cal_pi_k_tau();
    map<double,double> get_minimax_grid(const size_t &npoints);
    map<double,double> get_minimax_grid_tau(const size_t &npoints);
    //void cal_MP2_energy( map<double,double> &freq_grid );
   
    void cal_MP2_energy_pi_k_tau( map<double,double> &tau_grid ,std::map<double,map<Vector3_Order<double>,map<size_t,map<size_t,ComplexMatrix>>>> pi_k);

    //std::map<size_t,map<size_t,ComplexMatrix>> adj_atom_pair; 
    //void adj_atom_search();
*/ 
    const string minimax_grid_path="/home/shirong/LibRPA/minimax_grid";
    void cal_MP2_energy_pi_k( map<double,double> &freq_grid ); 
    void RPA_correlation_energy( map<double,double> &freq_grid );
    vector<double> read_cosine_trans_grid(const string &file_path);
    void dgeev(matrix &mat, vector<complex<double>> &egValue, matrix &vr, int info);
    int temp_Cs_count;  
    int temp_pi_count;
    int temp_chi0_k_count;
    double first_tau;
    double first_freq;
public:
    vector<ComplexMatrix> wfc_k;
    void chi0_main();
    void print_matrix( char* desc, const matrix &mat );
    void print_complex_matrix( char* desc, const ComplexMatrix &mat );
    void print_complex_real_matrix( char* desc, const ComplexMatrix &mat );
};

#endif 