#pragma once

#include <string>
#include <vector>

#include "../../src/core/meanfield.h"
#include "../../src/mpi/base_blacs.h"
#include "../../src/mpi/kpoint_blacs_parallel_context.h"

#include "reader_context.h"
namespace librpa::reader
{

enum class LegacyTextWfcOrder
{
    BasisSpinorBandSpin,
    SpinBasisBand
};

int read_eigenvector(ReaderContext &ctx, const std::string &dir_path);
int read_eigenvector(ReaderContext &ctx, const std::string &dir_path, librpa_int::MeanField &mf, bool use_spinor_wfc,
                     const std::vector<int> *iks_selected = nullptr);
int read_eigenvector(ReaderContext &ctx, const std::string &dir_path, librpa_int::MeanField &mf, bool use_spinor_wfc,
                     const std::vector<int> &source_to_target_ik,
                     const std::vector<int> *source_iks_selected,
                     LegacyTextWfcOrder text_order = LegacyTextWfcOrder::BasisSpinorBandSpin);

int read_eigenvector_kblacs_2d(ReaderContext &ctx,
    const std::string &dir_path, librpa_int::MeanField &mf, bool use_spinor_wfc,
    const librpa_int::KPointBlacsParallelContext &kblacs_ctxt,
    const librpa_int::ArrayDesc &desc_wfc,
    const std::vector<int> *source_to_target_ik = nullptr,
    LegacyTextWfcOrder text_order = LegacyTextWfcOrder::BasisSpinorBandSpin);

}
