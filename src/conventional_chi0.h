#pragma once

#include <vector>
#include <map>
#include "matrix.h"
#include "complexmatrix.h"
#include "vector3_order.h"
#include <memory>
#include "meanfield.h"
#include "timefreq.h"
#include "atoms.h"
#include "ri.h"
//#include "cal_periodic_chi0.h"
//#include "lib_main.h"
//#include "../src_global/global_variable.h"
class Conventional_Chi0
{
    //friend class Exx_Lcao;
//private:
    //const int ik=0;
    int range_all;
    double first_freq;
    //map<size_t,ComplexMatrix> C_tilde_map;  
    public:
    MeanField mf;
        const vector<Vector3_Order<double>> &klist;
        TFGrids tfg;
    Conventional_Chi0(const MeanField &mf_in,
             const vector<Vector3_Order<double>> &klist_in,
             unsigned n_tf_grids):
            mf(mf_in), klist(klist_in), tfg(n_tf_grids) {};
    ~Conventional_Chi0() {};
    void build(const Cs_LRI &Cs,
                   const vector<Vector3_Order<int>> &Rlist,
                   const Vector3_Order<int> &R_period,
                   const vector<atpair_t> &atpair_ABF,
                   const vector<Vector3_Order<double>> &qlist,
                   TFGrids::GRID_TYPES gt);
    map<size_t,ComplexMatrix> cal_ap_chi0_mat(const double &freq ,const size_t &Mu, const size_t &Nu,const vector<Vector3_Order<double>> &qlist);
    map<size_t,map<size_t,ComplexMatrix>> C_band_aux();
    complex<double> cal_cRPA_omega(map<size_t,map<size_t,ComplexMatrix>> &pi_mat, const double &freq, const vector<atom_t> &part_range);
    map<size_t,map<size_t,ComplexMatrix>> cal_pi_mat(map<size_t,map<size_t,ComplexMatrix>> &chi0_omega_mat, const double &freq);
    // map<size_t,ComplexMatrix> cal_ap_chi0_mat(const atpair_R_mat_t &Cs_IJR,
    //                                      const Vector3_Order<int> &R_period, 
    //                                      int spin_channel,
    //                                      atom_t mu, atom_t nu, double freq, vector<Vector3_Order<int>> qlist);
    // map<size_t,map<size_t,ComplexMatrix>> C_band_aux();

    map<size_t,map<size_t,map<size_t,map<size_t,ComplexMatrix>>>> aims_Cc(const Cs_LRI &Cs, const vector<Vector3_Order<double>> &qlist);   //Cc[is][ik][I][i_state](i,mu)
    void aims_C_tilde(const Cs_LRI &Cs, const vector<Vector3_Order<double>> &qlist); //aims_C_tilde[is][k][q][I](mn,mu)
    // ComplexMatrix construct_cc_plus_cc(const size_t &is, const size_t &I, const size_t &J);
    matrix construct_omega_band_mat(const size_t &is, const size_t ik, const size_t iq, const double &omega);
    
    // ComplexMatrix reshape_complexmat(const size_t n1, const size_t n2, const size_t n3, const ComplexMatrix &mat);
    // map<double,double> construct_gauss_grid(const int Npoints);
    // void print_matrix( char* desc, const matrix &mat );
    // void print_complex_matrix( char* desc, const ComplexMatrix &mat );
    // void print_complex_real_matrix( char* desc, const ComplexMatrix &mat );
    // vector<int> construct_part_range();
    
    // void ph_cRPA(const map<size_t,map<size_t,ComplexMatrix>> &C_tilde_map);
    // matrix construct_band_mat(const size_t &is);
    // matrix construct_delta_occ_mat(const size_t &is);
    // void dgeev(matrix &mat, vector<complex<double>> &egValue, matrix &vr, int info);
    
    // double get_band_gap();
    // map<double, double> get_minimax_grid(const size_t &npoints);
    // map<double,double> read_file_grid(const string &file_path, const char type);
    // //ComplexMatrix reshape_complexmat(const size_t n1, const size_t n2, const size_t n3, const ComplexMatrix &mat);
    // vector<pair<size_t,size_t>> get_atom_pair(const map<size_t,map<size_t,map<Vector3_Order<double>,std::shared_ptr<ComplexMatrix>>>> &m);
    size_t Tot_band_num;
    size_t occ_num[2];
    size_t unocc_num[2];
    size_t tot_ia[2];
    size_t tmp_count;
    map<size_t,map<size_t,map<size_t,map<size_t,ComplexMatrix>>>> C_tilde_map;
    // //size_t begin_occ;
    // // std::vector<size_t> occ_num;
    // // std::vector<size_t> unocc_num;
    // // std::vector<size_t> tot_ia; //occ_num*unocc_num
    // //std::map<size_t,size_t> occ_num;
    // //std::map<size_t,size_t> unocc_num;
    // //std::map<size_t,size_t> tot_ia; //occ_num*unocc_num
    // //std::vector<int> test_a;
    // vector<ComplexMatrix> wfc_k;
    
};

