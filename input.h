#ifndef INPUT_H
#define INPUT_H

#include "matrix.h"
#include <string>
#include "vector3.h"
#include "matrix3.h"
#include "complexmatrix.h"
#include<vector>
#include "vector3_order.h"
#include <map>
#include<memory>
extern const double PI;
extern const double TWO_PI;


extern int NBANDS;
extern int NLOCAL;
extern int NSPIN;
extern int n_kpoints;
extern int n_irk_points;
extern double efermi;
extern int natom;
extern int ncell;

extern matrix wg;
extern double **ekb;
extern int kv_nmp[3];
extern Matrix3 latvec;
extern Matrix3 G;
extern Vector3<double> *kvec_c;

extern map<int,int> atom_nw;
extern map<int,int> atom_mu;

extern map<Vector3_Order<double>,double>  irk_weight;
//extern vector<ComplexMatrix> wfc_k;	
extern map<size_t,map<size_t,map<Vector3_Order<int>,std::shared_ptr<matrix>>>> Cs;
extern map<size_t,map<size_t,map<Vector3_Order<double>,std::shared_ptr<ComplexMatrix>>>> Vq;
int atom_iw_loc2glo(const int &atom_index,const int & iw_lcoal);
int atom_mu_loc2glo(const int &atom_index,const int & mu_lcoal);
void READ_AIMS_BAND(const std::string &file_path);
void READ_AIMS_STRU(const std::string &file_path);
void READ_AIMS_EIGENVECTOR(const std::string &file_path,vector<ComplexMatrix> &wfc_k);
void READ_AIMS_Cs(const std::string &file_path);
void READ_AIMS_Vq(const std::string &file_path);
void handle_Cs_file(const std::string &file_path);
void handle_Vq_file(const std::string &file_path, map<Vector3_Order<double>,ComplexMatrix> &Vq_full);
void handle_KS_file(const std::string &file_path,vector<ComplexMatrix> &wfc_k);
#endif



