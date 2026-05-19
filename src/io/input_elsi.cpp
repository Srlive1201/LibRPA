#include "input_elsi.h"

#include <iostream>
#include <regex>

namespace librpa_int
{

// 读取 ELSI CSC 文件的函数
bool read_elsi_to_csc(const std::string& file_path, std::vector<int>& col_ptr,
                      std::vector<int>& row_idx, std::vector<std::complex<double>>& nnz_val,
                      std::vector<std::complex<double>>& nnz_val_cplx, int& n_basis,
                      const bool force_cplx)
{
    std::ifstream ifs(file_path, std::ios::binary);
    if (!ifs)
        throw LIBRPA_RUNTIME_ERROR("Cannot open file " + file_path);

    bool is_complex = true;

    // 读取文件内容
    ifs.seekg(0, std::ios::end);
    std::streampos size = ifs.tellg();
    ifs.seekg(0, std::ios::beg);

    std::vector<char> buffer(size);
    ifs.read(buffer.data(), size);
    ifs.close();

    // parse header
    int64_t header[16];
    std::memcpy(header, buffer.data(), 128);

    n_basis = header[3];
    int64_t nnz = header[5];

    // column
    int64_t* col_ptr_raw = reinterpret_cast<int64_t*>(buffer.data() + 128);
    col_ptr.assign(col_ptr_raw, col_ptr_raw + n_basis);
    col_ptr.push_back(nnz + 1);

    // row indices
    int32_t* row_idx_raw = reinterpret_cast<int32_t*>(buffer.data() + 128 + n_basis * 8);
    row_idx.assign(row_idx_raw, row_idx_raw + nnz);

    // non-zero values
    char* nnz_val_raw = buffer.data() + 128 + n_basis * 8 + nnz * 4;
    if (header[2] == 0)
    {
        double* nnz_ptr = reinterpret_cast<double*>(nnz_val_raw);
        if (force_cplx)
        {
            is_complex = true;
            nnz_val_cplx.resize(nnz);
            for (int64_t i = 0; i < nnz; ++i)
            {
                nnz_val_cplx[i] = std::complex<double>(nnz_ptr[i], 0.0);
            }
        }
        else
        {
            is_complex = false;
            nnz_val.assign(nnz_ptr, nnz_ptr + nnz);
        }
    }
    else
    {
        is_complex = true;
        auto *nnz_ptr_cplx = reinterpret_cast<std::complex<double> *>(nnz_val_raw);
        nnz_val_cplx.assign(nnz_ptr_cplx, nnz_ptr_cplx + nnz);
    }

    // Convert indices to 0-based
    for (int32_t& idx : row_idx) {
        idx -= 1;
    }
    for (int32_t& ptr : col_ptr) {
        ptr -= 1;
    }

    return is_complex;
}

static Matz load_matrix(const std::string& file_path, const MAJOR major = MAJOR::COL)
{
    std::vector<int> col_ptr;
    std::vector<int> row_idx;
    std::vector<double> nnz_val;
    std::vector<std::complex<double>> nnz_val_cplx;
    int n_basis;
    bool is_compelx;

    is_compelx = read_elsi_to_csc(file_path, col_ptr, row_idx, nnz_val, nnz_val_cplx, n_basis, true);
    return load_csc_to_matrix(n_basis, col_ptr, row_idx, nnz_val_cplx, major);
}

bool convert_csc(const std::string& filePath, std::map<std::string, Matz>& matrices, std::string& key) {
    // 定义正则表达式，用于匹配两种文件名格式
    std::regex filePattern(R"((\w+)_spin_(\d+)_kpt_(\d{6})(?:_freq_(\d+))?.csc|band_vxc_mat_spin_(\d+)_k_(\d{5}).csc)");
    std::smatch match;

    // 如果文件名符合正则表达式模式
    if (std::regex_search(filePath, match, filePattern)) {
        std::ostringstream oss;

        // 检查是否匹配到第二种格式
        if (!match[5].str().empty()) {
            // 第二种格式
            oss << "band_vxc_mat_spin_" << match[5] << "_k_" << match[6];
        } else {
            // 第一种格式
            oss << match[1] << "_spin_" << match[2] << "_kpt_" << std::setw(6) << std::setfill('0') << match[3];
            if (!match[4].str().empty()) {
                oss << "_freq_" << match[4];
            }
        }

        key = oss.str(); // 生成的key

        // 尝试加载矩阵，并处理可能的异常
        try {
            Matz matrix = load_matrix(filePath);
            matrices[key] = matrix;

            std::cout << "Matrix loaded and stored successfully under key: " << key << std::endl;
        } catch (const std::exception& e) {
            std::cerr << "Failed to load matrix from file: " << filePath << " Error: " << e.what() << std::endl;
            return false;
        }

        return true;
    } else {
        // 如果文件名不符合正则表达式模式，返回错误
        std::cerr << "无法解析文件名: " << filePath << std::endl;
        return false;
    }
}

}
