#include "vxc_io.h"

#include "../../src/io/fs.h"
#include "../../src/io/input_elsi.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace librpa_int
{
namespace qsgw
{
namespace
{

constexpr double hermitian_tolerance = 1.0e-10;

std::string trim(std::string value)
{
    value.erase(
        value.begin(),
        std::find_if(value.begin(), value.end(), [](const unsigned char ch) {
            return !std::isspace(ch);
        }));
    value.erase(
        std::find_if(value.rbegin(), value.rend(), [](const unsigned char ch) {
            return !std::isspace(ch);
        }).base(),
        value.end());
    return value;
}

std::string lowercase(std::string value)
{
    std::transform(value.begin(), value.end(), value.begin(),
                   [](const unsigned char ch) {
                       return static_cast<char>(std::tolower(ch));
                   });
    return value;
}

bool finite_complex(const cplxdb value)
{
    return std::isfinite(value.real()) && std::isfinite(value.imag());
}

void validate_hermitian_matrix(const Matz& matrix,
                               const int expected_dimension,
                               const std::string& label)
{
    if (matrix.nr() != expected_dimension ||
        matrix.nc() != expected_dimension)
    {
        throw std::invalid_argument(label + " has an invalid shape");
    }
    for (int row = 0; row < matrix.nr(); ++row)
    {
        for (int column = 0; column < matrix.nc(); ++column)
        {
            if (!finite_complex(matrix(row, column)))
            {
                throw std::invalid_argument(label + " contains non-finite data");
            }
            if (std::abs(matrix(row, column) -
                         std::conj(matrix(column, row))) >
                hermitian_tolerance)
            {
                throw std::invalid_argument(label + " is not Hermitian");
            }
        }
    }
}

std::vector<cplxdb> parse_complex_values(const std::string& line,
                                         const std::string& source_name)
{
    std::vector<cplxdb> result;
    std::size_t position = 0;
    while (true)
    {
        const std::size_t left = line.find('(', position);
        if (left == std::string::npos) break;
        const std::size_t comma = line.find(',', left + 1);
        const std::size_t right =
            line.find(')', comma == std::string::npos ? left + 1 : comma + 1);
        if (comma == std::string::npos || right == std::string::npos)
        {
            throw std::invalid_argument(
                "Malformed ABACUS complex value in " + source_name);
        }
        try
        {
            const double real = std::stod(
                trim(line.substr(left + 1, comma - left - 1)));
            const double imaginary = std::stod(
                trim(line.substr(comma + 1, right - comma - 1)));
            result.emplace_back(real, imaginary);
        }
        catch (const std::exception&)
        {
            throw std::invalid_argument(
                "Invalid ABACUS complex value in " + source_name);
        }
        position = right + 1;
    }
    return result;
}

Matz read_abacus_vxc_ha(std::istream& input,
                        const std::string& source_name)
{
    int rows = -1;
    int columns = -1;
    int current_row = -1;
    bool saw_row_marker = false;
    bool legacy_dimension_header = false;
    std::map<int, std::vector<cplxdb>> triangular_rows;
    std::vector<cplxdb> dense_values;
    std::string line;
    while (std::getline(input, line))
    {
        const std::string content = trim(line);
        if (content.empty()) continue;
        const std::string lowered = lowercase(content);
        if (lowered.rfind("# rows", 0) == 0)
        {
            if (legacy_dimension_header)
                throw std::invalid_argument(
                    "Mixed ABACUS Vxc matrix headers in " + source_name);
            std::istringstream fields(content);
            std::string comment_marker;
            std::string key;
            if (!(fields >> comment_marker >> key >> rows))
                throw std::invalid_argument(
                    "Malformed ABACUS Vxc row header in " + source_name);
            continue;
        }
        if (lowered.rfind("# columns", 0) == 0)
        {
            if (legacy_dimension_header)
                throw std::invalid_argument(
                    "Mixed ABACUS Vxc matrix headers in " + source_name);
            std::istringstream fields(content);
            std::string comment_marker;
            std::string key;
            if (!(fields >> comment_marker >> key >> columns))
                throw std::invalid_argument(
                    "Malformed ABACUS Vxc column header in " + source_name);
            continue;
        }
        if (content.front() == '#') continue;
        if (rows < 0 && columns < 0 && dense_values.empty() &&
            triangular_rows.empty() && content.find('(') == std::string::npos)
        {
            std::istringstream fields(content);
            int dimension = -1;
            if ((fields >> dimension) && dimension > 0 &&
                (fields >> std::ws).eof())
            {
                rows = dimension;
                columns = dimension;
                legacy_dimension_header = true;
                continue;
            }
        }
        if (lowered.rfind("row ", 0) == 0)
        {
            if (legacy_dimension_header)
                throw std::invalid_argument(
                    "Mixed ABACUS Vxc triangular layouts in " + source_name);
            std::istringstream fields(content);
            std::string label;
            int row_one_based = 0;
            if (!(fields >> label >> row_one_based) || row_one_based <= 0)
            {
                throw std::invalid_argument(
                    "Malformed ABACUS Vxc row marker in " + source_name);
            }
            current_row = row_one_based - 1;
            if (!triangular_rows.emplace(current_row,
                                         std::vector<cplxdb>{}).second)
            {
                throw std::invalid_argument(
                    "Duplicate ABACUS Vxc row marker in " + source_name);
            }
            saw_row_marker = true;
        }

        const std::vector<cplxdb> values =
            parse_complex_values(content, source_name);
        if (values.empty()) continue;
        if (saw_row_marker)
        {
            if (current_row < 0 || !dense_values.empty())
            {
                throw std::invalid_argument(
                    "Mixed ABACUS Vxc dense and triangular layouts in " +
                    source_name);
            }
            auto& row = triangular_rows.at(current_row);
            row.insert(row.end(), values.begin(), values.end());
        }
        else
        {
            dense_values.insert(dense_values.end(), values.begin(), values.end());
        }
    }

    if (rows <= 0 || columns <= 0 || rows != columns)
    {
        throw std::invalid_argument(
            "ABACUS Vxc matrix header is incomplete or non-square in " +
            source_name);
    }
    Matz result(rows, columns, MAJOR::ROW);
    if (legacy_dimension_header)
    {
        const std::size_t expected =
            static_cast<std::size_t>(rows) * (rows + 1) / 2;
        if (saw_row_marker || dense_values.size() != expected)
        {
            throw std::invalid_argument(
                "Legacy ABACUS triangular Vxc entry count mismatch in " +
                source_name);
        }
        std::size_t index = 0;
        for (int row = 0; row < rows; ++row)
        {
            for (int column = row; column < columns; ++column)
            {
                const cplxdb value = 0.5 * dense_values[index++];
                result(row, column) = value;
                if (row != column)
                    result(column, row) = std::conj(value);
            }
        }
    }
    else if (saw_row_marker)
    {
        if (!dense_values.empty() ||
            triangular_rows.size() != static_cast<std::size_t>(rows))
        {
            throw std::invalid_argument(
                "Incomplete ABACUS triangular Vxc matrix in " + source_name);
        }
        for (int row = 0; row < rows; ++row)
        {
            const auto row_it = triangular_rows.find(row);
            const int expected = columns - row;
            if (row_it == triangular_rows.end() ||
                static_cast<int>(row_it->second.size()) != expected)
            {
                throw std::invalid_argument(
                    "ABACUS triangular Vxc row length mismatch in " +
                    source_name);
            }
            for (int offset = 0; offset < expected; ++offset)
            {
                const int column = row + offset;
                const cplxdb value = 0.5 * row_it->second[offset];
                result(row, column) = value;
                if (row != column)
                    result(column, row) = std::conj(value);
            }
        }
    }
    else
    {
        if (dense_values.size() !=
            static_cast<std::size_t>(rows) * columns)
        {
            throw std::invalid_argument(
                "ABACUS dense Vxc entry count mismatch in " + source_name);
        }
        for (int row = 0; row < rows; ++row)
        {
            for (int column = 0; column < columns; ++column)
            {
                result(row, column) =
                    0.5 * dense_values[static_cast<std::size_t>(row) *
                                           columns + column];
            }
        }
    }
    validate_hermitian_matrix(result, rows, "ABACUS Vxc matrix");
    return result;
}

} // namespace

SpinKMatrixMap read_qsgw_vxc(const std::string& input_directory,
                            const std::string& prefix,
                            const MeanField& reference,
                            const bool aims_input,
                            const bool band_path,
                            const VxcBasis basis)
{
    if (!aims_input && reference.get_n_spinor() != 1)
        throw std::invalid_argument("ABACUS QSGW Vxc currently requires n_spinor=1");
    if (aims_input && basis != VxcBasis::State)
        throw std::invalid_argument("FHI-aims QSGW Vxc requires the KS-state basis");
    const std::string default_prefix = aims_input
        ? (band_path ? "band_vxc_mat" : "xc_matr")
        : (band_path ? "band_vxc" : "vxc");
    const std::string file_prefix = prefix.empty() ? default_prefix : prefix;
    const std::string path_prefix = is_absolute_path(file_prefix)
        ? file_prefix : join_path(input_directory, file_prefix);

    SpinKMatrixMap result;
    for (int spin = 0; spin < reference.get_n_spins(); ++spin)
        for (int kpoint = 0; kpoint < reference.get_n_kpoints(); ++kpoint)
        {
            std::ostringstream filename;
            filename << path_prefix;
            if (aims_input)
            {
                filename << "_spin_" << spin + 1
                         << (band_path ? "_k_" : "_kpt_")
                         << std::setw(band_path ? 5 : 6) << std::setfill('0')
                         << kpoint + 1 << ".csc";
            }
            else
            {
                filename << "k" << kpoint + 1;
                if (reference.get_n_spins() == 2) filename << "s" << spin + 1;
                filename << "_nao.txt";
            }
            std::string path = filename.str();
            // ABACUS Gamma-only exports omit the k-point suffix.
            if (!aims_input && reference.get_n_kpoints() == 1 &&
                !file_exists(path))
            {
                path = path_prefix;
                if (reference.get_n_spins() == 2) path += "s" + std::to_string(spin + 1);
                path += "_nao.txt";
            }
            require_readable_file(path);
            Matz matrix;
            if (aims_input)
                matrix = load_matrix_cplx(path, MAJOR::COL);
            else
            {
                std::ifstream input(path);
                matrix = read_abacus_vxc_ha(input, path);
            }
            result[spin][kpoint] = prepare_vxc_in_fixed_state_basis(
                matrix, basis, reference, spin, kpoint);
        }
    return result;
}

} // namespace qsgw
} // namespace librpa_int
