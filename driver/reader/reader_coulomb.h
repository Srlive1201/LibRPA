#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include "../../src/core/ri.h"

#include "reader_context.h"
namespace librpa::reader
{

bool check_coulomb_file_binary(const std::string &file_path);

int detect_coulomb_reader_version(const std::string &dir_path,
                                  const std::string &vq_fprefix);

size_t read_Vq_full(ReaderContext &ctx, const std::string &dir_path, const std::string &vq_fprefix,
                    bool is_cut_coulomb, int reader_version = 0,
                    bool use_shrink_basis = false);

size_t read_Vq_row(ReaderContext &ctx, const std::string &dir_path, const std::string &vq_fprefix,
                   double threshold, const std::vector<librpa_int::atpair_t> &local_atpair,
                   bool is_cut_coulomb, int reader_version = 0,
                   bool use_shrink_basis = false);

}
