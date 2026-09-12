//! Functions to parse LibRPA driver input files
/*
 */
#ifndef READ_DATA_H
#define READ_DATA_H

#include <string>
#include <vector>

#include "../src/core/ri.h"
#include "../src/math/matrix.h"
#include "../src/math/vector3_order.h"
#include "librpa.hpp"
#include "reader_eigenvec.h"
#include "reader_lri.h"
#include "reader_coulomb.h"

// TODO: remove this include and internal datatypes in signature.
// Data objects of internal types should be accessed in the implementation.
#include "../src/core/meanfield.h"
#include "../src/core/dielecmodel.h"

using std::string;
using librpa_int::matrix;
using librpa_int::atpair_t;
using librpa_int::Vector3_Order;
using librpa_int::atpair_R_mat_t;
using librpa_int::MeanField;
using librpa_int::ComplexMatrix;
using librpa_int::velocity_matrix_t;

/*!
 * @brief Read occupation numbers and eigenvalues of SCF calculation
 */
void read_scf_occ_eigenvalues(const string &file_path);
void read_scf_occ_eigenvalues(const string &file_path, MeanField &mf, bool use_spinor_wfc);
void read_scf_occ_eigenvalues(const string &file_path, MeanField &mf, bool use_spinor_wfc,
                              const std::vector<int> &source_ik_for_target,
                              int source_n_kpoints);

/*!
 * @brief Read exchange-correlation potential
 *
 * @param[in] file_path   data file path
 * @return    status code, 0 for succesful read, 1 with error
 */
int read_vxc(const string &file_path, std::vector<matrix> &vxc);

// high-level reader for RI coefficients and bare Coulomb interactions
void read_ri(const string &dir_path, librpa::ParallelRouting &routing);

void read_velocity(const string &file_path, const MeanField &mf, velocity_matrix_t &velocity);
void read_velocity(const string &file_path, const MeanField &mf, velocity_matrix_t &velocity,
                   const std::vector<int> &source_to_target_ik, int source_n_kpoints);
void read_velocity_abacus(const MeanField &mf, const std::string &dir_path,
                          const std::string &file_prefix, velocity_matrix_t &velocity);
void read_velocity_aims(const MeanField &mf, const std::string &file_path,
                        velocity_matrix_t &velocity);
void read_velocity_aims(const MeanField &mf, const std::string &file_path,
                        const std::string &file_prefix, velocity_matrix_t &velocity);
void read_headwing_input(const string &dir_path, bool need_wing);

void read_ri_shrink(const string &dir_path);

size_t read_shrink_sinvS(const string &dir_path, const string &vq_fprefix,
                         std::map<Vector3_Order<double>, ComplexMatrix> &sinvS);

void read_stru(const std::string &file_path);

void read_bz_sampling(const std::string &file_path);
void read_bz_sampling_from_stru(const std::string &file_path);

void read_basis_wfc_aux(const std::string &input_dir, const std::string &fn_basis,
                        const std::string &fn_basis_wfc, const std::string &fn_basis_aux);

void read_dielec_func(const string &file_path, std::vector<double> &omegas,
                      std::vector<double> &dielec_func_imagfreq);

void erase_Cs_from_local_atp(atpair_R_mat_t &Cs, std::vector<atpair_t> &local_atpair);

void read_band_kpath_info(const string &file_path);

void read_band_meanfield_data(const string &dir_path);

std::vector<matrix> read_vxc_band(const string &dir_path, int n_states, int n_spin,
                                  int n_kpoints_band);

void read_elsi_csc(const std::string &file_path, bool save_row_major, std::vector<double> &mat,
                   int &n_basis, bool &is_real);
#endif  // !READ_DATA_H
