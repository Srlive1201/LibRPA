#pragma once

#include <istream>
#include <map>
#include <string>
#include <vector>

#include "../../src/math/complexmatrix.h"
#include "../../src/math/vector3_order.h"

namespace librpa::reader
{

struct NaoOverlapRealSpace
{
    int dimension = 0;
    std::map<librpa_int::Vector3_Order<int>, librpa_int::ComplexMatrix> blocks;
};

//! Read ABACUS verbose S(R) CSR or compact get_s CSR with an optional STEP preamble.
NaoOverlapRealSpace read_abacus_nao_overlap_csr(std::istream& input,
                                                const std::string& source_name = "<stream>");

//! Open and read an ABACUS S(R) CSR file.
NaoOverlapRealSpace read_abacus_nao_overlap_csr(const std::string& path);

//! Fold S(R) using the ABACUS convention S(k)=sum_R exp(+2*pi*i*k.R) S(R).
std::vector<librpa_int::ComplexMatrix> fourier_nao_overlap(
    const NaoOverlapRealSpace& overlap,
    const std::vector<librpa_int::Vector3_Order<double>>& kfrac_list);

//! Validate/extract S from an interleaved ABACUS I_spinor tensor S export.
//! A spatial export is also accepted. No spin averaging conceals a discrepancy.
librpa_int::ComplexMatrix spatial_nao_overlap(const librpa_int::ComplexMatrix& overlap, int n_aos,
                                              int n_spinor);

}  // namespace librpa::reader
