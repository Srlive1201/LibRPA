#include "convert_csc.h"

#include <iostream>
#include <regex>
#include <filesystem>


// 读取 ELSI CSC 文件的函数
static void read_elsi_to_csc(const std::string& filePath, std::vector<int>& col_ptr, std::vector<int>& row_idx, std::vector<std::complex<double>>& nnz_val, int& n_basis) {
    std::ifstream inputFile(filePath, std::ios::binary);
    if (!inputFile) {
        std::cerr << "无法打开文件: " << filePath << std::endl;
        return;
    }

    // 读取文件内容
    inputFile.seekg(0, std::ios::end);
    std::streampos size = inputFile.tellg();
    inputFile.seekg(0, std::ios::beg);

    std::vector<char> buffer(size);
    inputFile.read(buffer.data(), size);
    inputFile.close();

    // 解析 header
    int64_t header[16];
    std::memcpy(header, buffer.data(), 128);

    n_basis = header[3];
    int64_t nnz = header[5];

    // 获取列指针
    int64_t* col_ptr_raw = reinterpret_cast<int64_t*>(buffer.data() + 128);
    col_ptr.assign(col_ptr_raw, col_ptr_raw + n_basis);
    col_ptr.push_back(nnz + 1);

    // 获取行索引
    int32_t* row_idx_raw = reinterpret_cast<int32_t*>(buffer.data() + 128 + n_basis * 8);
    row_idx.assign(row_idx_raw, row_idx_raw + nnz);

    // 获取非零值
    char* nnz_val_raw = buffer.data() + 128 + n_basis * 8 + nnz * 4;
    if (header[2] == 0) {
        // 实数情况
        double* nnz_val_double = reinterpret_cast<double*>(nnz_val_raw);
        for (int64_t i = 0; i < nnz; ++i) {
            nnz_val.push_back(std::complex<double>(nnz_val_double[i], 0.0));
        }
    } else {
        // 复数情况
        double* nnz_val_double = reinterpret_cast<double*>(nnz_val_raw);
        for (int64_t i = 0; i < nnz; ++i) {
            nnz_val.push_back(std::complex<double>(nnz_val_double[2 * i], nnz_val_double[2 * i + 1]));
        }
    }

    // 更改索引
    for (int32_t& idx : row_idx) {
        idx -= 1;
    }
    for (int32_t& ptr : col_ptr) {
        ptr -= 1;
    }
}


static Matz loadMatrix(const std::string& filePath) {
    std::vector<int> col_ptr;
    std::vector<int> row_idx;
    std::vector<std::complex<double>> nnz_val;
    int n_basis;

    read_elsi_to_csc(filePath, col_ptr, row_idx, nnz_val, n_basis);

    Matz matrix(n_basis, n_basis, MAJOR::COL);
    for (int col = 0; col < n_basis; ++col) {
        for (int idx = col_ptr[col]; idx < col_ptr[col + 1]; ++idx) {
            int row = row_idx[idx];
            matrix(row, col) = nnz_val[idx];
        }
    }
    

    return matrix;
}

// bool convert_csc(const std::string& filePath, std::map<std::string, Matz>& matrices, std::string& key)
// {
//     // 定义正则表达式，用于匹配文件名格式
//     std::regex filePattern(R"((\w+)_spin_(\d+)_kpt_(\d{6})(?:_freq_(\d+))?.csc)");
//     std::smatch match;

//     // 如果文件名符合正则表达式模式
//     if (std::regex_search(filePath, match, filePattern)) {
//         std::string fileType = match[1];
//         int ispin = std::stoi(match[2]) - 1;  // 自旋索引从1基转换为0基
//         int ikpt = std::stoi(match[3]) - 1;   // k点从1基转换为0基
//         int freq = match[4].matched ? std::stoi(match[4]) : -1;

//         std::cout << "File Type: " << fileType << std::endl;
//         std::cout << "Spin: " << ispin << std::endl;
//         std::cout << "K-point: " << ikpt << std::endl;
//         if (freq != -1) {
//             std::cout << "Frequency: " << freq << std::endl;
//         }

//         std::cout << "尝试打开文件: " << filePath << std::endl;

//         // 加载矩阵
//         Matz matrix = loadMatrix(filePath);

//         // 使用ostringstream格式化key，包括6位的k点
//         std::ostringstream oss;
//         oss << fileType << "_spin_" << (ispin + 1) << "_kpt_" << std::setw(6) << std::setfill('0') << (ikpt + 1);
//         if (freq != -1) {
//             oss << "_freq_" << freq;
//         }

//         key = oss.str();  // 生成的key
//         // 将矩阵存储在 map 中
//         matrices[key] = matrix;

//         std::cout << "Matrix loaded and stored successfully under key: " << key << std::endl;
//         return true;
//     }

//     // 如果文件名不符合正则表达式模式，返回错误
//     std::cerr << "无法解析文件名: " << filePath << std::endl;
//     return false;
// }

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
            Matz matrix = loadMatrix(filePath);
            
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