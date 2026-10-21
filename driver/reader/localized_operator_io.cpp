#include "localized_operator_io.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <locale>
#include <set>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <vector>

#include "../../src/utils/constants.h"

namespace librpa::reader
{
namespace
{
using namespace librpa_int;

std::string line(std::istream& input)
{
    std::string result;
    while (std::getline(input, result))
    {
        const auto comment = result.find('#');
        if (comment != std::string::npos) result.erase(comment);
        if (result.find_first_not_of(" \t\r") != std::string::npos) return result;
    }
    throw std::runtime_error("Incomplete QSGW AO CSR input");
}

int last_integer(const std::string& text)
{
    std::istringstream input(text);
    std::string token;
    while (input >> token)
    {
    }
    std::size_t end = 0;
    const int value = std::stoi(token, &end);
    if (end != token.size()) throw std::runtime_error("Invalid QSGW AO CSR header");
    return value;
}

enum class CsrFormat
{
    Compact,
    Legacy,
    Modern,
};

std::tuple<int, int, CsrFormat> header(std::istream& input)
{
    std::string first = line(input);
    if (first.find("--- Ionic Step") != std::string::npos)
    {
        int n = -1, nr = -1;
        std::string text;
        while (std::getline(input, text))
        {
            if (text.find("--- Ionic Step") != std::string::npos)
                throw std::runtime_error("Multiple ionic steps in QSGW AO CSR input");
            if (text.find("number of localized basis") != std::string::npos)
            {
                if (!(std::istringstream(text) >> n))
                    throw std::runtime_error("Invalid QSGW AO CSR basis count");
            }
            if (text.find("number of Bravais lattice vector R") != std::string::npos)
            {
                if (!(std::istringstream(text) >> nr))
                    throw std::runtime_error("Invalid QSGW AO CSR R count");
            }
            if (text.find("CSR Format") != std::string::npos)
            {
                if (n <= 0 || nr < 0) throw std::runtime_error("Invalid QSGW AO CSR dimensions");
                return {n, nr, CsrFormat::Modern};
            }
        }
        throw std::runtime_error("Missing CSR Format section in QSGW AO input");
    }
    const bool native = first.find("STEP:") != std::string::npos ||
                        first.find("Matrix Dimension") != std::string::npos;
    if (first.find("STEP:") != std::string::npos) first = line(input);
    const int n = last_integer(first);
    const int nr = last_integer(line(input));
    if (n <= 0 || nr < 0) throw std::runtime_error("Invalid QSGW AO CSR dimensions");
    return {n, nr, native ? CsrFormat::Legacy : CsrFormat::Compact};
}

void skip_comment_lines(std::istream& input)
{
    input >> std::ws;
    while (input.peek() == '#')
    {
        std::string ignored;
        std::getline(input, ignored);
        input >> std::ws;
    }
}

void csr_marker(std::istream& input, const std::string& expected)
{
    std::string text;
    while (std::getline(input, text))
    {
        const auto first = text.find_first_not_of(" \t\r");
        if (first == std::string::npos) continue;
        const auto last = text.find_last_not_of(" \t\r");
        if (text.substr(first, last - first + 1) != expected)
            throw std::runtime_error("Expected " + expected + " in QSGW AO CSR input");
        return;
    }
    throw std::runtime_error("Missing " + expected + " in QSGW AO CSR input");
}

Matz matrix(std::istream& input, const int n, const int nnz, const bool empty_has_pointers,
            const bool upper_only = true, const CsrFormat format = CsrFormat::Legacy)
{
    if (nnz < 0 || static_cast<long long>(nnz) > static_cast<long long>(n) * n)
        throw std::runtime_error("Invalid QSGW AO CSR nonzero count");
    Matz result(n, n, MAJOR::ROW);
    std::fill(result.ptr(), result.ptr() + result.size(), cplxdb{});
    if (nnz == 0 && !empty_has_pointers && format != CsrFormat::Modern) return result;
    std::vector<cplxdb> values(nnz);
    std::vector<int> columns(nnz), pointers(n + 1);
    if (format == CsrFormat::Modern) csr_marker(input, "# CSR values");
    for (auto& x : values)
        if (!(input >> x) || !std::isfinite(x.real()) || !std::isfinite(x.imag()))
            throw std::runtime_error("Invalid QSGW AO CSR value");
    if (format == CsrFormat::Modern) csr_marker(input, "# CSR column indices");
    for (auto& x : columns)
        if (!(input >> x) || x < 0 || x >= n)
            throw std::runtime_error("Invalid QSGW AO CSR column");
    if (format == CsrFormat::Modern) csr_marker(input, "# CSR row pointers");
    for (auto& x : pointers)
        if (!(input >> x) || x < 0 || x > nnz)
            throw std::runtime_error("Invalid QSGW AO CSR row pointer");
    if (pointers.front() != 0 || pointers.back() != nnz ||
        !std::is_sorted(pointers.begin(), pointers.end()))
        throw std::runtime_error("Invalid QSGW AO CSR row range");
    for (int row = 0; row < n; ++row)
    {
        int previous = -1;
        for (int j = pointers[row]; j < pointers[row + 1]; ++j)
        {
            if (columns[j] <= previous)
                throw std::runtime_error("QSGW AO CSR columns are not unique and ordered");
            previous = columns[j];
            if (!upper_only || columns[j] >= row) result(row, columns[j]) = values[j];
        }
    }
    return result;
}

std::set<Vector3_Order<int>> csr_support(const librpa_int::LocalizedMatrix& blocks, const int n,
                                         const double scale)
{
    if (n <= 0 || n > std::numeric_limits<int>::max() / n)
        throw std::runtime_error("Invalid QSGW Hamiltonian CSR dimension");
    std::set<Vector3_Order<int>> support;
    for (const auto& [r, block] : blocks)
    {
        if (block.nr() != n || block.nc() != n)
            throw std::runtime_error("Invalid QSGW Hamiltonian CSR shape");
        for (std::size_t j = 0; j < block.size(); ++j)
        {
            const auto value = block.ptr()[j] * scale;
            if (!std::isfinite(value.real()) || !std::isfinite(value.imag()))
                throw std::runtime_error("Non-finite QSGW Hamiltonian CSR value");
        }
        support.insert(r);
    }
    return support;
}

void write_full_hamiltonian_csr(std::ostream& output, const librpa_int::LocalizedMatrix& blocks,
                                const std::set<Vector3_Order<int>>& support, const int n,
                                const double scale, const int iteration, const int precision,
                                const bool modern = false, const int spin_directions = 1,
                                const int spin_index = 1)
{
    if (modern)
    {
        output << " --- Ionic Step " << iteration << " ---\n"
               << " # print H matrix in real space H(R)\n"
               << ' ' << spin_directions << " # number of spin directions\n " << spin_index
               << " # spin index\n " << n << " # number of localized basis\n " << support.size()
               << " # number of Bravais lattice vector R\n\n # CSR Format\n";
    }
    else
    {
        output << "STEP: " << iteration << "\nMatrix Dimension of H(R): " << n
               << "\nMatrix number of H(R): " << support.size() << '\n';
    }
    output << std::setprecision(precision) << std::scientific;
    for (const auto& r : support)
    {
        const auto found = blocks.find(r);
        // Explicit full CSR also represents all-zero padded blocks, avoiding
        // the consumers' different zero-nnz payload conventions.
        output << r.x << ' ' << r.y << ' ' << r.z << ' ' << n * n << '\n';
        if (modern) output << " # CSR values\n";
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                output << (found == blocks.end() ? cplxdb{} : found->second(i, j) * scale) << ' ';
        output << '\n';
        if (modern) output << " # CSR column indices\n";
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) output << j << ' ';
        output << '\n';
        if (modern) output << " # CSR row pointers\n";
        for (int i = 0; i <= n; ++i) output << i * n << ' ';
        output << '\n';
    }
}
}  // namespace

librpa_int::LocalizedAoOperators read_headwing_ao_operators(const std::string& overlap_file,
                                                            const std::string& position_file,
                                                            const double position_scale_bohr)
{
    if (!(position_scale_bohr > 0.0) || !std::isfinite(position_scale_bohr))
        throw std::runtime_error("Invalid head/wing position-unit conversion");
    std::ifstream s(overlap_file), r(position_file);
    if (!s || !r) throw std::runtime_error("Cannot open QSGW head/wing S(R) or r(R) CSR input");
    const auto [n, ns, format_s] = header(s);
    const auto [nr, nr_blocks, format_r] = header(r);
    if (n != nr) throw std::runtime_error("QSGW head/wing S(R) and r(R) dimensions differ");
    librpa_int::LocalizedAoOperators result;
    result.dimension = n;
    for (int i = 0; i < ns; ++i)
    {
        int x, y, z, nnz;
        skip_comment_lines(s);
        if (!(s >> x >> y >> z >> nnz)) throw std::runtime_error("Incomplete QSGW overlap R block");
        if (!result.overlap
                 .emplace(Vector3_Order<int>(x, y, z), matrix(s, n, nnz, false, true, format_s))
                 .second)
            throw std::runtime_error("Duplicate QSGW overlap R block");
    }
    for (int i = 0; i < nr_blocks; ++i)
    {
        int x, y, z;
        skip_comment_lines(r);
        if (!(r >> x >> y >> z)) throw std::runtime_error("Incomplete QSGW position R block");
        for (int a = 0; a < 3; ++a)
        {
            int nnz;
            if (!(r >> nnz)) throw std::runtime_error("Incomplete QSGW position nonzero count");
            if (!result.position[a]
                     .emplace(Vector3_Order<int>(x, y, z),
                              matrix(r, n, nnz, format_r == CsrFormat::Compact, true, format_r))
                     .second)
                throw std::runtime_error("Duplicate QSGW position R block");
        }
    }
    std::string trailing;
    skip_comment_lines(s);
    skip_comment_lines(r);
    if (s >> trailing || r >> trailing) throw std::runtime_error("Trailing QSGW AO CSR data");
    for (auto& component : result.position)
        for (auto& [cell, block] : component)
            for (std::size_t j = 0; j < block.size(); ++j) block.ptr()[j] *= position_scale_bohr;
    return result;
}

librpa_int::LocalizedMatrix read_headwing_hamiltonian_csr(const std::string& filename,
                                                          const double energy_scale_ha)
{
    if (!(energy_scale_ha > 0.0) || !std::isfinite(energy_scale_ha))
        throw std::runtime_error("Invalid head/wing Hamiltonian-unit conversion");
    std::ifstream input(filename);
    if (!input) throw std::runtime_error("Cannot open head/wing H(R) CSR input: " + filename);
    const auto [n, nr, format] = header(input);
    librpa_int::LocalizedMatrix result;
    for (int i = 0; i < nr; ++i)
    {
        int x, y, z, nnz;
        skip_comment_lines(input);
        if (!(input >> x >> y >> z >> nnz))
            throw std::runtime_error("Incomplete head/wing Hamiltonian R block");
        auto block = matrix(input, n, nnz, false, false, format);
        for (std::size_t j = 0; j < block.size(); ++j) block.ptr()[j] *= energy_scale_ha;
        if (!result.emplace(Vector3_Order<int>(x, y, z), std::move(block)).second)
            throw std::runtime_error("Duplicate head/wing Hamiltonian R block");
    }
    std::string trailing;
    skip_comment_lines(input);
    if (input >> trailing) throw std::runtime_error("Trailing head/wing Hamiltonian CSR data");
    return result;
}
void write_headwing_hamiltonian_csr(const std::string& filename,
                                    const librpa_int::LocalizedMatrix& blocks)
{
    if (blocks.empty()) throw std::runtime_error("Empty QSGW head/wing Hamiltonian");
    const int n = blocks.begin()->second.nr();
    const auto support = csr_support(blocks, n, librpa_int::HA2EV);
    std::ofstream output(filename);
    if (!output) throw std::runtime_error("Cannot write QSGW head/wing Hamiltonian: " + filename);
    write_full_hamiltonian_csr(output, blocks, support, n, librpa_int::HA2EV, 0, 17);
    output.flush();
    if (!output) throw std::runtime_error("Cannot flush QSGW head/wing Hamiltonian: " + filename);
}

void write_qsgw_hamiltonian_csr(const std::string& filename,
                                const librpa_int::LocalizedMatrix& blocks,
                                const librpa_int::LocalizedAoOperators& operators,
                                const int iteration, const int spin_directions,
                                const int spin_index)
{
    if (blocks.empty()) throw std::runtime_error("Empty QSGW Hamiltonian for export");
    if (iteration < 0) throw std::runtime_error("Invalid QSGW Hamiltonian export iteration");
    if ((spin_directions != 1 && spin_directions != 2) || spin_index < 1 ||
        spin_index > spin_directions)
        throw std::runtime_error("Invalid QSGW Hamiltonian export spin channel");
    const int n = blocks.begin()->second.nr();
    if (operators.dimension != n || operators.overlap.empty())
        throw std::runtime_error("QSGW Hamiltonian export AO operators differ from H(R)");
    auto support = csr_support(blocks, n, 2.0);
    const auto add_support = [&](const librpa_int::LocalizedMatrix& localized)
    {
        const auto extra = csr_support(localized, n, 1.0);
        support.insert(extra.begin(), extra.end());
    };
    add_support(operators.overlap);
    for (const auto& component : operators.position) add_support(component);

    const std::string temporary = filename + ".tmp";
    try
    {
        std::ofstream output(temporary);
        if (!output) throw std::runtime_error("Cannot write QSGW Hamiltonian: " + temporary);
        output.imbue(std::locale::classic());
        // Scientific precision 16 gives 17 significant decimal digits.
        write_full_hamiltonian_csr(output, blocks, support, n, 2.0, iteration, 16, true,
                                   spin_directions, spin_index);
        output.flush();
        if (!output) throw std::runtime_error("Cannot flush QSGW Hamiltonian: " + temporary);
        output.close();
        if (!output) throw std::runtime_error("Cannot close QSGW Hamiltonian: " + temporary);
        std::filesystem::rename(temporary, filename);
    }
    catch (...)
    {
        std::error_code ignored;
        std::filesystem::remove(temporary, ignored);
        throw;
    }
}
}  // namespace librpa::reader
