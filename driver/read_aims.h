//! Functions to parse files generated by FHI-aims
/*
 */
#ifndef READ_AIMS_H
#define READ_AIMS_H

#include <string>
#include <map>
#include <set>
#include "meanfield.h"
#include "ri.h"
#include "vector3_order.h"
using std::string;
using std::map;

void READ_AIMS_BAND(const string &file_path, MeanField &mf);
void READ_AIMS_EIGENVECTOR(const string &dir_path, MeanField &mf);
void handle_KS_file(const string &file_path, MeanField &mf);
size_t READ_AIMS_Cs(const string &dir_path, double threshold,const vector<atpair_t> &local_atpair);
size_t READ_AIMS_Cs_evenly_distribute(const string &dir_path, double threshold, int myid, int nprocs);
size_t READ_Vq_Full(const string &dir_path, const string &vq_fprefix, double threshold, atpair_k_cplx_mat_t &coulomb);
size_t READ_Vq_Row(const string &dir_path, const string &vq_fprefix, double threshold, atpair_k_cplx_mat_t &coulomb, const vector<atpair_t> &local_atpair);
void READ_AIMS_STRU(const int& n_kpoints, const std::string &file_path);

void read_dielec_func(const string &file_path, std::vector<double> &omegas, std::vector<double> &dielec_func_imagfreq);
size_t handle_Cs_file(const std::string &file_path, double threshold,const vector<atpair_t> &local_atpair);

std::vector<size_t> handle_Cs_file_dry(const string &file_path, double threshold);

size_t handle_Cs_file_by_ids(const string &file_path, double threshold, const vector<size_t> &ids);

void handle_Vq_full_file(const std::string &file_path, double threshold, map<Vector3_Order<double>, ComplexMatrix> &Vq_full);
void handle_Vq_row_file(const std::string &file_path, double threshold, atpair_k_cplx_mat_t &coulomb, const vector<atpair_t> &local_atpair);
void erase_Cs_from_local_atp(atpair_R_mat_t &Cs, vector<atpair_t> &local_atpair);

void read_aims(MeanField &mf);
size_t get_natom_ncell_from_first_Cs_file(const string file_path = "Cs_data_0.txt");

#endif // !READ_AIMS_H