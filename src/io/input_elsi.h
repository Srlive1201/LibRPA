#pragma once

#include <map>
#include <string>

#include "../math/matrix_m.h"

namespace librpa_int
{

//! Read ELSI CSC file into standard CSC format
bool read_elsi_to_csc(const std::string& filePath, std::vector<int>& col_ptr,
                      std::vector<int>& row_idx,
                      std::vector<double>& nnz_val,
                      std::vector<std::complex<double>>& nnz_val_cplx,
                      int& n_basis, const bool force_cplx = true);

template <typename T>
matrix_m<T> load_csc_to_matrix(const int n_basis, const std::vector<int>& col_ptr, const std::vector<int>& row_idx,
                               const std::vector<T>& nnz_val, const MAJOR major = MAJOR::COL)
{
    matrix_m<T> mat(n_basis, n_basis, major);
    for (int col = 0; col < n_basis; ++col)
    {
        for (int idx = col_ptr[col]; idx < col_ptr[col + 1]; ++idx)
        {
            int row = row_idx[idx];
            mat(row, col) = nnz_val[idx];
        }
    }

    return mat;
}

bool convert_csc(const std::string& filePath, std::map<std::string, Matz>& matrices, std::string& key);

}  // namespace librpa_int
