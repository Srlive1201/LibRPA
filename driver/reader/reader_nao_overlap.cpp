#include "reader_nao_overlap.h"

#include <cmath>
#include <complex>
#include <fstream>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace librpa::reader
{

librpa_int::ComplexMatrix spatial_nao_overlap(const librpa_int::ComplexMatrix& overlap, int n_aos,
                                              int n_spinor)
{
    if (n_aos <= 0 || (n_spinor != 1 && n_spinor != 2) || overlap.nr != overlap.nc)
        throw std::runtime_error("invalid cRPA overlap dimensions");
    if (overlap.nr == n_aos) return overlap;
    if (n_spinor != 2 || overlap.nr != 2 * n_aos)
        throw std::runtime_error("cRPA overlap and loaded spatial AO dimensions differ");
    librpa_int::ComplexMatrix spatial(n_aos, n_aos);
    for (int i = 0; i < n_aos; ++i)
        for (int j = 0; j < n_aos; ++j)
        {
            spatial(i, j) = overlap(2 * i, 2 * j);
            const double tolerance = 1e-12 * (1 + std::abs(spatial(i, j)));
            if (std::abs(overlap(2 * i + 1, 2 * j + 1) - spatial(i, j)) > tolerance ||
                std::abs(overlap(2 * i, 2 * j + 1)) > tolerance ||
                std::abs(overlap(2 * i + 1, 2 * j)) > tolerance)
                throw std::runtime_error("cRPA exported overlap is not spin independent");
        }
    return spatial;
}

namespace
{

struct LineCursor
{
    std::vector<std::string> lines;
    std::size_t index = 0;
    std::string source;
};

std::string trim(const std::string& text)
{
    const std::size_t begin = text.find_first_not_of(" \t\r\n");
    if (begin == std::string::npos) return {};
    const std::size_t end = text.find_last_not_of(" \t\r\n");
    return text.substr(begin, end - begin + 1);
}

bool is_comment(const std::string& line)
{
    const std::string value = trim(line);
    return !value.empty() && value.front() == '#';
}

[[noreturn]] void fail(const LineCursor& cursor, const std::string& message)
{
    std::ostringstream out;
    out << cursor.source << ':' << (cursor.index + 1) << ": " << message;
    throw std::runtime_error(out.str());
}

std::string next_payload(LineCursor& cursor)
{
    while (cursor.index < cursor.lines.size())
    {
        const std::string line = trim(cursor.lines[cursor.index++]);
        if (!line.empty() && !is_comment(line)) return line;
    }
    fail(cursor, "unexpected end of file");
}

long long parse_integer_token(const std::string& token, const LineCursor& cursor,
                              const std::string& context)
{
    std::size_t consumed = 0;
    long long value = 0;
    try
    {
        value = std::stoll(token, &consumed);
    }
    catch (const std::exception&)
    {
        fail(cursor, "invalid integer in " + context + ": " + token);
    }
    if (consumed != token.size()) fail(cursor, "invalid integer in " + context + ": " + token);
    return value;
}

int parse_labeled_integer(const std::string& line, const std::string& label,
                          const LineCursor& cursor)
{
    const std::size_t label_pos = line.find(label);
    if (label_pos == std::string::npos) fail(cursor, "missing header label: " + label);

    std::string payload = trim(line.substr(0, label_pos));
    if (!payload.empty() && payload.back() == '#')
        payload = trim(payload.substr(0, payload.size() - 1));
    if (payload.empty())
    {
        payload = trim(line.substr(label_pos + label.size()));
        if (!payload.empty() && payload.front() == ':') payload = trim(payload.substr(1));
    }

    std::istringstream parser(payload);
    std::string token;
    std::string last;
    while (parser >> token) last = token;
    if (last.empty()) fail(cursor, "missing integer for header label: " + label);
    const long long value = parse_integer_token(last, cursor, label);
    if (value < 0 || value > std::numeric_limits<int>::max())
        fail(cursor, "integer out of range for " + label);
    return static_cast<int>(value);
}

std::size_t find_unique_line(const LineCursor& cursor, const std::string& marker)
{
    std::size_t found = cursor.lines.size();
    for (std::size_t i = 0; i < cursor.lines.size(); ++i)
    {
        if (cursor.lines[i].find(marker) == std::string::npos) continue;
        if (found != cursor.lines.size())
        {
            LineCursor copy = cursor;
            copy.index = i;
            fail(copy, "duplicate header marker: " + marker);
        }
        found = i;
    }
    if (found == cursor.lines.size())
    {
        LineCursor copy = cursor;
        copy.index = 0;
        fail(copy, "missing header marker: " + marker);
    }
    return found;
}

bool consume_optional_comment_marker(LineCursor& cursor, const std::string& marker)
{
    while (cursor.index < cursor.lines.size() && trim(cursor.lines[cursor.index]).empty())
        ++cursor.index;
    if (cursor.index >= cursor.lines.size()) fail(cursor, "unexpected end before " + marker);
    const std::string line = trim(cursor.lines[cursor.index]);
    if (!is_comment(line)) return false;
    ++cursor.index;
    if (line.find(marker) == std::string::npos)
        fail(cursor, "expected comment marker " + marker + ", got: " + line);
    return true;
}

void expect_comment_marker(LineCursor& cursor, const std::string& marker)
{
    if (!consume_optional_comment_marker(cursor, marker))
        fail(cursor, "missing comment marker: " + marker);
}

std::vector<std::string> read_tokens(LineCursor& cursor, std::size_t count,
                                     const std::string& context)
{
    std::vector<std::string> result;
    result.reserve(count);
    while (result.size() < count)
    {
        if (cursor.index >= cursor.lines.size())
            fail(cursor, "unexpected end while reading " + context);
        const std::string line = trim(cursor.lines[cursor.index++]);
        if (line.empty()) continue;
        if (is_comment(line)) fail(cursor, "too few values in " + context);

        std::istringstream parser(line);
        std::string token;
        while (parser >> token)
        {
            if (result.size() == count) fail(cursor, "too many values in " + context);
            result.push_back(token);
        }
    }
    return result;
}

std::complex<double> parse_complex_token(const std::string& token, const LineCursor& cursor)
{
    std::complex<double> value;
    std::istringstream parser(token);
    parser >> value;
    if (!parser) fail(cursor, "invalid finite overlap value: " + token);
    parser >> std::ws;
    if (parser.peek() != std::char_traits<char>::eof() || !std::isfinite(value.real()) ||
        !std::isfinite(value.imag()))
        fail(cursor, "invalid finite overlap value: " + token);
    return value;
}

bool parse_header(LineCursor& cursor, int& dimension, int& number_of_r)
{
    std::size_t first = 0;
    while (first < cursor.lines.size() && trim(cursor.lines[first]).empty()) ++first;
    if (first == cursor.lines.size()) fail(cursor, "empty overlap file");

    std::size_t compact_header = first;
    std::string first_line = trim(cursor.lines[compact_header]);
    if (first_line.rfind("STEP:", 0) == 0)
    {
        LineCursor step_cursor = cursor;
        step_cursor.index = compact_header;
        const std::string step_token = trim(first_line.substr(5));
        const long long step = parse_integer_token(step_token, step_cursor, "STEP");
        if (step < 0) fail(step_cursor, "STEP must be non-negative");
        ++compact_header;
        while (compact_header < cursor.lines.size() && trim(cursor.lines[compact_header]).empty())
            ++compact_header;
        if (compact_header == cursor.lines.size())
            fail(step_cursor, "unexpected end after STEP header");
        first_line = trim(cursor.lines[compact_header]);
    }
    if (first_line.find("Matrix Dimension") != std::string::npos)
    {
        cursor.index = compact_header + 1;
        // Native ABACUS uses its generic H-container CSR header for H and S.
        // The S(k) metric/spin checks below still validate the overlap itself.
        dimension = parse_labeled_integer(first_line, "Matrix Dimension", cursor);
        const std::string count_line = next_payload(cursor);
        number_of_r = parse_labeled_integer(count_line, "Matrix number", cursor);
        return true;
    }

    const std::size_t title = find_unique_line(cursor, "print S matrix in real space S(R)");
    const std::size_t nspin = find_unique_line(cursor, "number of spin directions");
    const std::size_t spin = find_unique_line(cursor, "spin index");
    const std::size_t basis = find_unique_line(cursor, "number of localized basis");
    const std::size_t rcount = find_unique_line(cursor, "number of Bravais lattice vector R");
    const std::size_t csr = find_unique_line(cursor, "CSR Format");
    if (!(first < title && title < nspin && nspin < spin && spin < basis && basis < rcount &&
          rcount < csr))
        fail(cursor, "verbose overlap header fields are out of order");

    if (parse_labeled_integer(cursor.lines[nspin], "# number of spin directions", cursor) != 1)
        fail(cursor, "S(R) must contain exactly one spin-independent overlap channel");
    if (parse_labeled_integer(cursor.lines[spin], "# spin index", cursor) != 1)
        fail(cursor, "S(R) spin index must be one");
    dimension = parse_labeled_integer(cursor.lines[basis], "# number of localized basis", cursor);
    number_of_r =
        parse_labeled_integer(cursor.lines[rcount], "# number of Bravais lattice vector R", cursor);
    cursor.index = csr + 1;
    return false;
}

std::vector<int> parse_integer_block(LineCursor& cursor, std::size_t count,
                                     const std::string& context)
{
    const auto tokens = read_tokens(cursor, count, context);
    std::vector<int> values;
    values.reserve(count);
    for (const auto& token : tokens)
    {
        const long long value = parse_integer_token(token, cursor, context);
        if (value < std::numeric_limits<int>::min() || value > std::numeric_limits<int>::max())
            fail(cursor, "integer out of range in " + context);
        values.push_back(static_cast<int>(value));
    }
    return values;
}

}  // namespace

NaoOverlapRealSpace read_abacus_nao_overlap_csr(std::istream& input, const std::string& source_name)
{
    LineCursor cursor;
    cursor.source = source_name;
    std::string line;
    while (std::getline(input, line)) cursor.lines.push_back(line);
    if (!input.eof()) throw std::runtime_error(source_name + ": failed while reading overlap file");

    int number_of_r = 0;
    NaoOverlapRealSpace result;
    const bool compact_format = parse_header(cursor, result.dimension, number_of_r);
    if (result.dimension <= 0) fail(cursor, "overlap dimension must be positive");
    if (number_of_r <= 0) fail(cursor, "overlap file must contain at least one R block");

    for (int block_index = 0; block_index < number_of_r; ++block_index)
    {
        const std::string block_header = next_payload(cursor);
        std::istringstream header_parser(block_header);
        std::string rx_token, ry_token, rz_token, nnz_token, trailing;
        if (!(header_parser >> rx_token >> ry_token >> rz_token >> nnz_token) ||
            (header_parser >> trailing))
            fail(cursor, "R-block header must contain exactly four integers");

        const long long rx_long = parse_integer_token(rx_token, cursor, "R x");
        const long long ry_long = parse_integer_token(ry_token, cursor, "R y");
        const long long rz_long = parse_integer_token(rz_token, cursor, "R z");
        const long long nnz_long = parse_integer_token(nnz_token, cursor, "R block nnz");
        if (rx_long < std::numeric_limits<int>::min() ||
            rx_long > std::numeric_limits<int>::max() ||
            ry_long < std::numeric_limits<int>::min() ||
            ry_long > std::numeric_limits<int>::max() ||
            rz_long < std::numeric_limits<int>::min() || rz_long > std::numeric_limits<int>::max())
            fail(cursor, "R coordinate is out of int range");
        const long long maximum_nnz = static_cast<long long>(result.dimension) * result.dimension;
        if (nnz_long < 0 || nnz_long > maximum_nnz)
            fail(cursor, "R-block nnz is outside [0, dimension^2]");

        const librpa_int::Vector3_Order<int> r(static_cast<int>(rx_long), static_cast<int>(ry_long),
                                               static_cast<int>(rz_long));
        if (result.blocks.find(r) != result.blocks.end()) fail(cursor, "duplicate R block");

        if (compact_format && nnz_long == 0)
        {
            result.blocks.emplace(r, librpa_int::ComplexMatrix(result.dimension, result.dimension));
            continue;
        }

        const std::size_t nnz = static_cast<std::size_t>(nnz_long);
        const bool has_section_markers = consume_optional_comment_marker(cursor, "CSR values");
        const auto value_tokens = read_tokens(cursor, nnz, "CSR values");
        if (has_section_markers) expect_comment_marker(cursor, "CSR column indices");
        const auto columns = parse_integer_block(cursor, nnz, "CSR column indices");
        if (has_section_markers) expect_comment_marker(cursor, "CSR row pointers");
        const auto rows = parse_integer_block(
            cursor, static_cast<std::size_t>(result.dimension + 1), "CSR row pointers");

        if (rows.front() != 0 || rows.back() != static_cast<int>(nnz))
            fail(cursor, "CSR row pointers must start at zero and end at nnz");
        for (int row = 0; row < result.dimension; ++row)
        {
            if (rows[row] > rows[row + 1] || rows[row] < 0 || rows[row + 1] > static_cast<int>(nnz))
                fail(cursor, "CSR row pointers are not monotone and bounded");
            int previous_column = -1;
            for (int offset = rows[row]; offset < rows[row + 1]; ++offset)
            {
                if (columns[offset] < 0 || columns[offset] >= result.dimension)
                    fail(cursor, "CSR column index is out of range");
                if (columns[offset] <= previous_column)
                    fail(cursor, "CSR columns must be strictly increasing within each row");
                previous_column = columns[offset];
            }
        }

        librpa_int::ComplexMatrix dense(result.dimension, result.dimension);
        for (int row = 0; row < result.dimension; ++row)
        {
            for (int offset = rows[row]; offset < rows[row + 1]; ++offset)
                dense(row, columns[offset]) = parse_complex_token(value_tokens[offset], cursor);
        }
        result.blocks.emplace(r, std::move(dense));
    }

    while (cursor.index < cursor.lines.size())
    {
        const std::string remainder = trim(cursor.lines[cursor.index++]);
        if (!remainder.empty() && !is_comment(remainder))
            fail(cursor, "unexpected trailing payload after declared R blocks");
    }
    if (static_cast<int>(result.blocks.size()) != number_of_r)
        fail(cursor, "parsed R-block count differs from header");
    return result;
}

NaoOverlapRealSpace read_abacus_nao_overlap_csr(const std::string& path)
{
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open ABACUS overlap file: " + path);
    return read_abacus_nao_overlap_csr(input, path);
}

std::vector<librpa_int::ComplexMatrix> fourier_nao_overlap(
    const NaoOverlapRealSpace& overlap,
    const std::vector<librpa_int::Vector3_Order<double>>& kfrac_list)
{
    if (overlap.dimension <= 0 || overlap.blocks.empty())
        throw std::invalid_argument("cannot Fourier transform an empty NAO overlap");

    constexpr double two_pi = 6.283185307179586476925286766559;
    std::vector<librpa_int::ComplexMatrix> result;
    result.reserve(kfrac_list.size());
    for (const auto& k : kfrac_list)
    {
        if (!std::isfinite(k.x) || !std::isfinite(k.y) || !std::isfinite(k.z))
            throw std::invalid_argument("NAO overlap k point must be finite");
        librpa_int::ComplexMatrix sk(overlap.dimension, overlap.dimension);
        for (const auto& item : overlap.blocks)
        {
            const auto& r = item.first;
            const double argument = two_pi * (k.x * r.x + k.y * r.y + k.z * r.z);
            const std::complex<double> phase(std::cos(argument), std::sin(argument));
            const auto& sr = item.second;
            if (sr.nr != overlap.dimension || sr.nc != overlap.dimension)
                throw std::invalid_argument("NAO overlap R block has inconsistent dimensions");
            for (int i = 0; i < overlap.dimension; ++i)
                for (int j = 0; j < overlap.dimension; ++j) sk(i, j) += phase * sr(i, j);
        }
        result.push_back(std::move(sk));
    }
    return result;
}

}  // namespace librpa::reader
