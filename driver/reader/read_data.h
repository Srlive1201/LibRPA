#pragma once
#include <string>
#include <vector>

#include "../../src/core/dielecmodel.h"
#include "../../src/core/meanfield.h"
#include "reader_context.h"
namespace librpa::reader
{
using librpa_int::MeanField;
using librpa_int::velocity_matrix_t;
using std::string;
void read_scf_occ_eigenvalues(const string &file_path, MeanField &mf, bool use_spinor_wfc);
void read_scf_occ_eigenvalues(const string &file_path, MeanField &mf, bool use_spinor_wfc,
                              const std::vector<int> &source_ik_for_target, int source_n_kpoints);
void read_velocity(ReaderContext &ctx, const string &file_path, const MeanField &mf,
                   velocity_matrix_t &velocity);
void read_velocity(ReaderContext &ctx, const string &file_path, const MeanField &mf,
                   velocity_matrix_t &velocity, const std::vector<int> &source_to_target_ik,
                   int source_n_kpoints);
void read_velocity_abacus(ReaderContext &ctx, const MeanField &mf, const std::string &dir_path,
                          const std::string &file_prefix, velocity_matrix_t &velocity);
void read_velocity_aims(ReaderContext &ctx, const MeanField &mf, const std::string &file_path,
                        const std::string &file_prefix, velocity_matrix_t &velocity);
void read_bz_sampling(ReaderContext &ctx, const std::string &file_path);
void read_bz_sampling_from_stru(ReaderContext &ctx, const std::string &file_path);
void read_basis_wfc_aux(ReaderContext &ctx, const std::string &input_dir,
                        const std::string &fn_basis, const std::string &fn_basis_wfc,
                        const std::string &fn_basis_aux);
void read_band_kpath_info(ReaderContext &ctx, const string &file_path);
void read_band_eigenvalues(ReaderContext &ctx, const string &dir_path);
void read_band_meanfield_data(ReaderContext &ctx, const string &dir_path);
}  // namespace librpa::reader
