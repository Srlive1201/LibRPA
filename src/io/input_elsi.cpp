#include "input_elsi.h"

#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <regex>

#include "fs.h"

namespace librpa_int
{

namespace
{

std::string get_file_name(const std::string& file_path)
{
    const auto pos = file_path.find_last_of("/\\");
    if (pos == std::string::npos) return file_path;
    return file_path.substr(pos + 1);
}

bool parse_regular_csc_key(const std::string& file_name, std::string& key)
{
    static const std::regex pattern(
        R"(^([[:alnum:]_]+)_spin_([0-9]+)_kpt_([0-9]{6})(?:_freq_([0-9]+))?\.csc$)");

    std::smatch match;
    if (!std::regex_match(file_name, match, pattern)) return false;

    key = match[1].str() + "_spin_" + match[2].str() + "_kpt_" + match[3].str();
    if (match[4].matched) key += "_freq_" + match[4].str();
    return true;
}

bool parse_band_vxc_csc_key(const std::string& file_name, std::string& key)
{
    static const std::regex pattern(R"(^band_vxc_mat_spin_([0-9]+)_k_([0-9]{5})\.csc$)");

    std::smatch match;
    if (!std::regex_match(file_name, match, pattern)) return false;

    key = "band_vxc_mat_spin_" + match[1].str() + "_k_" + match[2].str();
    return true;
}

bool parse_csc_key(const std::string& file_path, std::string& key)
{
    const auto file_name = get_file_name(file_path);
    return parse_band_vxc_csc_key(file_name, key) || parse_regular_csc_key(file_name, key);
}

}  // namespace

// Read ELSI CSC file into standard CSC format
bool read_elsi_to_csc(const std::string& file_path, std::vector<int>& col_ptr,
                      std::vector<int>& row_idx, std::vector<double>& nnz_val,
                      std::vector<std::complex<double>>& nnz_val_cplx, int& n_basis,
                      const bool force_cplx)
{
    require_readable_file(file_path);
    std::ifstream ifs(file_path, std::ios::binary);
    if (!ifs) throw LIBRPA_RUNTIME_ERROR("Cannot open file " + file_path);

    // Read file content
    ifs.seekg(0, std::ios::end);
    const std::streamoff size = ifs.tellg();
    int64_t header[16];
    if (size < static_cast<std::streamoff>(sizeof(header)) ||
        static_cast<uint64_t>(size) > std::numeric_limits<std::size_t>::max())
        throw LIBRPA_RUNTIME_ERROR("Invalid ELSI CSC file size in " + file_path);
    ifs.seekg(0, std::ios::beg);

    std::vector<char> buffer(static_cast<std::size_t>(size));
    if (!ifs.read(buffer.data(), size))
        throw LIBRPA_RUNTIME_ERROR("Failed to read ELSI CSC file " + file_path);
    ifs.close();

    std::memcpy(header, buffer.data(), sizeof(header));
    if (header[2] != 0 && header[2] != 1)
        throw LIBRPA_RUNTIME_ERROR("Unsupported ELSI CSC value type in " + file_path);
    if (header[3] <= 0 || header[3] > std::numeric_limits<int>::max() || header[5] < 0 ||
        header[5] > std::numeric_limits<int>::max())
        throw LIBRPA_RUNTIME_ERROR("Invalid or unsupported ELSI CSC dimensions in " + file_path);

    const uint64_t column_bytes = static_cast<uint64_t>(header[3]) * sizeof(int64_t);
    const uint64_t row_bytes = static_cast<uint64_t>(header[5]) * sizeof(int32_t);
    const uint64_t value_bytes = static_cast<uint64_t>(header[5]) *
                                 (header[2] == 0 ? sizeof(double) : sizeof(std::complex<double>));
    if (sizeof(header) + column_bytes + row_bytes + value_bytes > buffer.size())
        throw LIBRPA_RUNTIME_ERROR("Truncated ELSI CSC matrix data in " + file_path);

    n_basis = static_cast<int>(header[3]);
    const int nnz = static_cast<int>(header[5]);

    // Validate file indices before converting them to zero-based indices.
    col_ptr.resize(static_cast<std::size_t>(n_basis) + 1);
    for (int col = 0; col < n_basis; ++col)
    {
        int64_t ptr;
        std::memcpy(&ptr,
                    buffer.data() + sizeof(header) + static_cast<std::size_t>(col) * sizeof(ptr),
                    sizeof(ptr));
        if (ptr < 1 || ptr > header[5] + 1 || (col == 0 && ptr != 1) ||
            (col > 0 && ptr - 1 < col_ptr[col - 1]))
            throw LIBRPA_RUNTIME_ERROR("Invalid ELSI CSC column pointer in " + file_path);
        col_ptr[col] = static_cast<int>(ptr - 1);
    }
    col_ptr[n_basis] = nnz;

    row_idx.resize(nnz);
    for (int i = 0; i < nnz; ++i)
    {
        int32_t row;
        std::memcpy(&row,
                    buffer.data() + sizeof(header) + column_bytes +
                        static_cast<std::size_t>(i) * sizeof(row),
                    sizeof(row));
        if (row < 1 || row > n_basis)
            throw LIBRPA_RUNTIME_ERROR("Invalid ELSI CSC row index in " + file_path);
        row_idx[i] = row - 1;
    }

    // non-zero values
    // Values can start at an unaligned offset when nnz is odd.
    const char* nnz_val_raw = buffer.data() + sizeof(header) + column_bytes + row_bytes;
    const bool is_complex = header[2] == 1 || force_cplx;
    if (header[2] == 0)
    {
        if (force_cplx)
        {
            nnz_val_cplx.resize(nnz);
            for (int i = 0; i < nnz; ++i)
            {
                double value;
                std::memcpy(&value, nnz_val_raw + static_cast<std::size_t>(i) * sizeof(value),
                            sizeof(value));
                nnz_val_cplx[i] = std::complex<double>(value, 0.0);
            }
        }
        else
        {
            nnz_val.resize(nnz);
            if (nnz > 0) std::memcpy(nnz_val.data(), nnz_val_raw, value_bytes);
        }
    }
    else
    {
        nnz_val_cplx.resize(nnz);
        if (nnz > 0) std::memcpy(nnz_val_cplx.data(), nnz_val_raw, value_bytes);
    }

    return is_complex;
}

Matz load_matrix_cplx(const std::string& file_path, const MAJOR major)
{
    std::vector<int> col_ptr;
    std::vector<int> row_idx;
    std::vector<double> nnz_val;
    std::vector<std::complex<double>> nnz_val_cplx;
    int n_basis;
    bool is_cplx;

    is_cplx = read_elsi_to_csc(file_path, col_ptr, row_idx, nnz_val, nnz_val_cplx, n_basis, true);
    if (is_cplx)
        return load_csc_to_matrix(n_basis, col_ptr, row_idx, nnz_val_cplx, major);
    return load_csc_to_matrix(n_basis, col_ptr, row_idx, nnz_val, major).to_complex();
}

Matd load_matrix_real(const std::string& file_path, const MAJOR major)
{
    std::vector<int> col_ptr;
    std::vector<int> row_idx;
    std::vector<double> nnz_val;
    std::vector<std::complex<double>> nnz_val_cplx;
    int n_basis;
    bool is_cplx;

    is_cplx = read_elsi_to_csc(file_path, col_ptr, row_idx, nnz_val, nnz_val_cplx, n_basis, true);
    if (is_cplx)
        return load_csc_to_matrix(n_basis, col_ptr, row_idx, nnz_val_cplx, major).get_real();
    return load_csc_to_matrix(n_basis, col_ptr, row_idx, nnz_val, major);
}

bool convert_csc(const std::string& filePath, std::map<std::string, Matz>& matrices,
                 std::string& key, const MAJOR major)
{
    if (!parse_csc_key(filePath, key))
    {
        std::cerr << "Failed to parse CSC file name: " << filePath << std::endl;
        return false;
    }

    try
    {
        matrices[key] = load_matrix_cplx(filePath, major);
        std::cout << "Matrix loaded and stored successfully under key: " << key << std::endl;
    }
    catch (const std::exception& e)
    {
        std::cerr << "Failed to load matrix from file: " << filePath << " Error: " << e.what()
                  << std::endl;
        return false;
    }

    return true;
}

}
