#include "localized_operator_io.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <vector>

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
    throw std::runtime_error("Incomplete Head/wing AO CSR input");
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
    if (end != token.size()) throw std::runtime_error("Invalid Head/wing AO CSR header");
    return value;
}

std::tuple<int, int, bool> header(std::istream& input)
{
    std::string first = line(input);
    const bool native = first.find("STEP:") != std::string::npos ||
                        first.find("Matrix Dimension") != std::string::npos;
    if (first.find("STEP:") != std::string::npos) first = line(input);
    const int n = last_integer(first);
    const int nr = last_integer(line(input));
    if (n <= 0 || nr < 0) throw std::runtime_error("Invalid Head/wing AO CSR dimensions");
    return {n, nr, native};
}

Matz matrix(std::istream& input, const int n, const int nnz, const bool empty_has_pointers,
            const bool upper_only = true)
{
    if (nnz < 0 || static_cast<long long>(nnz) > static_cast<long long>(n) * n)
        throw std::runtime_error("Invalid Head/wing AO CSR nonzero count");
    Matz result(n, n, MAJOR::ROW);
    std::fill(result.ptr(), result.ptr() + result.size(), cplxdb{});
    if (nnz == 0 && !empty_has_pointers) return result;
    std::vector<cplxdb> values(nnz);
    std::vector<int> columns(nnz), pointers(n + 1);
    for (auto& x : values)
        if (!(input >> x) || !std::isfinite(x.real()) || !std::isfinite(x.imag()))
            throw std::runtime_error("Invalid Head/wing AO CSR value");
    for (auto& x : columns)
        if (!(input >> x) || x < 0 || x >= n)
            throw std::runtime_error("Invalid Head/wing AO CSR column");
    for (auto& x : pointers)
        if (!(input >> x) || x < 0 || x > nnz)
            throw std::runtime_error("Invalid Head/wing AO CSR row pointer");
    if (pointers.front() != 0 || pointers.back() != nnz ||
        !std::is_sorted(pointers.begin(), pointers.end()))
        throw std::runtime_error("Invalid Head/wing AO CSR row range");
    for (int row = 0; row < n; ++row)
    {
        int previous = -1;
        for (int j = pointers[row]; j < pointers[row + 1]; ++j)
        {
            if (columns[j] <= previous)
                throw std::runtime_error("Head/wing AO CSR columns are not unique and ordered");
            previous = columns[j];
            if (!upper_only || columns[j] >= row) result(row, columns[j]) = values[j];
        }
    }
    return result;
}
}  // namespace

librpa_int::LocalizedAoOperators read_headwing_ao_operators(const std::string& overlap_file,
                                                            const std::string& position_file,
                                                            const double position_scale_bohr)
{
    if (!(position_scale_bohr > 0.0) || !std::isfinite(position_scale_bohr))
        throw std::runtime_error("Invalid head/wing position-unit conversion");
    std::ifstream s(overlap_file), r(position_file);
    if (!s || !r)
        throw std::runtime_error("Cannot open Head/wing head/wing S(R) or r(R) CSR input");
    const auto [n, ns, native_s] = header(s);
    const auto [nr, nr_blocks, native_r] = header(r);
    if (n != nr) throw std::runtime_error("Head/wing head/wing S(R) and r(R) dimensions differ");
    librpa_int::LocalizedAoOperators result;
    result.dimension = n;
    for (int i = 0; i < ns; ++i)
    {
        int x, y, z, nnz;
        if (!(s >> x >> y >> z >> nnz))
            throw std::runtime_error("Incomplete Head/wing overlap R block");
        if (!result.overlap.emplace(Vector3_Order<int>(x, y, z), matrix(s, n, nnz, false)).second)
            throw std::runtime_error("Duplicate Head/wing overlap R block");
    }
    for (int i = 0; i < nr_blocks; ++i)
    {
        int x, y, z;
        if (!(r >> x >> y >> z)) throw std::runtime_error("Incomplete Head/wing position R block");
        for (int a = 0; a < 3; ++a)
        {
            int nnz;
            if (!(r >> nnz))
                throw std::runtime_error("Incomplete Head/wing position nonzero count");
            if (!result.position[a]
                     .emplace(Vector3_Order<int>(x, y, z), matrix(r, n, nnz, !native_r))
                     .second)
                throw std::runtime_error("Duplicate Head/wing position R block");
        }
    }
    std::string trailing;
    if (s >> trailing || r >> trailing) throw std::runtime_error("Trailing Head/wing AO CSR data");
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
    const auto [n, nr, native] = header(input);
    librpa_int::LocalizedMatrix result;
    for (int i = 0; i < nr; ++i)
    {
        int x, y, z, nnz;
        if (!(input >> x >> y >> z >> nnz))
            throw std::runtime_error("Incomplete head/wing Hamiltonian R block");
        auto block = matrix(input, n, nnz, false, false);
        for (std::size_t j = 0; j < block.size(); ++j) block.ptr()[j] *= energy_scale_ha;
        if (!result.emplace(Vector3_Order<int>(x, y, z), std::move(block)).second)
            throw std::runtime_error("Duplicate head/wing Hamiltonian R block");
    }
    std::string trailing;
    if (input >> trailing) throw std::runtime_error("Trailing head/wing Hamiltonian CSR data");
    return result;
}
}  // namespace librpa::reader
