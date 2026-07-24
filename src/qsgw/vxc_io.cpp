#include "vxc_io.h"
#include "sha256.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <sstream>
#include <stdexcept>

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

void validate_reference_wfc(const MeanField& reference,
                            const int spin,
                            const int kpoint)
{
    if (!reference.initialized() || spin < 0 ||
        spin >= reference.get_n_spins() || kpoint < 0 ||
        kpoint >= reference.get_n_kpoints())
    {
        throw std::invalid_argument(
            "QSGW Vxc projection received an invalid mean-field index");
    }
    for (int spinor = 0; spinor < reference.get_n_spinor(); ++spinor)
    {
        const ComplexMatrix* wfc = reference.find_wfc(spin, spinor, kpoint);
        if (wfc == nullptr || wfc->nr != reference.get_n_bands() ||
            wfc->nc != reference.get_n_aos())
        {
            throw std::invalid_argument(
                "QSGW Vxc projection reference wavefunction is incomplete");
        }
        for (int index = 0; index < wfc->size; ++index)
        {
            if (!finite_complex(wfc->c[index]))
            {
                throw std::invalid_argument(
                    "QSGW Vxc projection wavefunction contains non-finite data");
            }
        }
    }
}

double periodic_distance(const double lhs, const double rhs)
{
    const double difference = lhs - rhs;
    return std::abs(difference - std::round(difference));
}

} // namespace

VxcManifest VxcManifest::parse(std::istream& input,
                               const std::string& source_name)
{
    VxcManifest result;
    std::map<std::string, std::string> metadata;
    bool saw_magic = false;
    bool saw_table_header = false;
    std::string line;
    while (std::getline(input, line))
    {
        const std::string content = trim(line);
        if (content.empty()) continue;
        if (!saw_magic)
        {
            if (content != "# librpa-qsgw-vxc-manifest-v2")
            {
                throw std::invalid_argument(
                    "Invalid QSGW Vxc manifest header in " + source_name);
            }
            saw_magic = true;
            continue;
        }
        if (!saw_table_header)
        {
            std::istringstream fields(content);
            std::string first;
            fields >> first;
            if (lowercase(first) == "spin")
            {
                std::vector<std::string> header;
                header.push_back(first);
                std::string field;
                while (fields >> field) header.push_back(field);
                const std::vector<std::string> expected{
                    "spin", "k_index", "kx", "ky", "kz", "rows",
                    "columns", "sha256", "file"};
                if (header.size() != expected.size())
                {
                    throw std::invalid_argument(
                        "Invalid QSGW Vxc manifest table header in " +
                        source_name);
                }
                for (std::size_t index = 0; index < expected.size(); ++index)
                {
                    if (lowercase(header[index]) != expected[index])
                    {
                        throw std::invalid_argument(
                            "Invalid QSGW Vxc manifest table header in " +
                            source_name);
                    }
                }
                saw_table_header = true;
                continue;
            }

            std::string value;
            fields >> value;
            std::string extra;
            if (first.empty() || value.empty() || fields >> extra)
            {
                throw std::invalid_argument(
                    "Malformed QSGW Vxc manifest metadata in " + source_name);
            }
            first = lowercase(first);
            if (!metadata.emplace(first, value).second)
            {
                throw std::invalid_argument(
                    "Duplicate QSGW Vxc manifest metadata in " + source_name);
            }
            continue;
        }

        VxcManifestEntry entry;
        std::istringstream fields(content);
        int spin_one_based = 0;
        int kpoint_one_based = 0;
        if (!(fields >> spin_one_based >> kpoint_one_based >> entry.kpoint.x >>
              entry.kpoint.y >> entry.kpoint.z >> entry.rows >>
              entry.columns >> entry.sha256 >> entry.file))
        {
            throw std::invalid_argument(
                "Malformed QSGW Vxc manifest entry in " + source_name);
        }
        std::string extra;
        if (fields >> extra || spin_one_based <= 0 || kpoint_one_based <= 0 ||
            entry.rows <= 0 || entry.columns <= 0 ||
            !is_sha256_hex(entry.sha256) || entry.file.empty() ||
            std::filesystem::path(entry.file).is_absolute() ||
            entry.file.find('\\') != std::string::npos ||
            !std::isfinite(entry.kpoint.x) ||
            !std::isfinite(entry.kpoint.y) ||
            !std::isfinite(entry.kpoint.z))
        {
            throw std::invalid_argument(
                "Invalid QSGW Vxc manifest entry in " + source_name);
        }
        for (const auto& component : std::filesystem::path(entry.file))
        {
            if (component == "..")
            {
                throw std::invalid_argument(
                    "QSGW Vxc manifest file escapes its base directory in " +
                    source_name);
            }
        }
        entry.spin = spin_one_based - 1;
        entry.k_index = kpoint_one_based - 1;
        if (!result.entries_
                 .emplace(std::make_pair(entry.spin, entry.k_index), entry)
                 .second)
        {
            throw std::invalid_argument(
                "Duplicate spin/k entry in QSGW Vxc manifest " + source_name);
        }
    }

    if (!saw_magic || !saw_table_header || result.entries_.empty())
    {
        throw std::invalid_argument(
            "Incomplete QSGW Vxc manifest " + source_name);
    }
    for (const char* key : {"kind", "producer", "units", "basis", "gauge"})
    {
        if (metadata.count(key) == 0)
        {
            throw std::invalid_argument(
                "Missing QSGW Vxc manifest metadata in " + source_name);
        }
    }
    if (metadata.size() != 5)
    {
        throw std::invalid_argument(
            "Unknown QSGW Vxc manifest metadata in " + source_name);
    }

    const std::string kind = lowercase(metadata.at("kind"));
    if (kind == "scf")
        result.kind_ = VxcDatasetKind::ScfGrid;
    else if (kind == "band")
        result.kind_ = VxcDatasetKind::BandPath;
    else
        throw std::invalid_argument(
            "Unsupported QSGW Vxc manifest kind in " + source_name);

    result.producer_ = lowercase(metadata.at("producer"));
    const std::string units = lowercase(metadata.at("units"));
    const std::string basis = lowercase(metadata.at("basis"));
    const std::string gauge = lowercase(metadata.at("gauge"));
    if (units == "ha" || units == "hartree")
        result.units_ = VxcUnits::Hartree;
    else if (units == "ry" || units == "rydberg")
        result.units_ = VxcUnits::Rydberg;
    else
        throw std::invalid_argument(
            "Unsupported QSGW Vxc units in " + source_name);
    if (basis == "nao")
        result.basis_ = VxcBasis::Nao;
    else if (basis == "state")
        result.basis_ = VxcBasis::State;
    else
        throw std::invalid_argument(
            "Unsupported QSGW Vxc basis in " + source_name);
    if (gauge == "ao_bloch")
        result.gauge_ = VxcGauge::AoBloch;
    else if (gauge == "mf0_state")
        result.gauge_ = VxcGauge::Mf0State;
    else
        throw std::invalid_argument(
            "Unsupported QSGW Vxc gauge in " + source_name);

    const bool valid_abacus_state =
        result.producer_ == "abacus" &&
        result.units_ == VxcUnits::Rydberg &&
        result.basis_ == VxcBasis::State &&
        result.gauge_ == VxcGauge::Mf0State;
    const bool valid_abacus_nao =
        result.producer_ == "abacus" &&
        result.units_ == VxcUnits::Rydberg &&
        result.basis_ == VxcBasis::Nao &&
        result.gauge_ == VxcGauge::AoBloch;
    const bool valid_aims =
        result.producer_ == "fhi-aims" &&
        result.units_ == VxcUnits::Hartree &&
        result.basis_ == VxcBasis::State &&
        result.gauge_ == VxcGauge::Mf0State;
    if (!valid_abacus_state && !valid_abacus_nao && !valid_aims)
    {
        throw std::invalid_argument(
            "Incompatible QSGW Vxc producer, units, basis, or gauge in " +
            source_name);
    }
    return result;
}

void VxcManifest::validate(
    const VxcDatasetKind expected_kind,
    const int spin_count,
    const std::vector<Vector3_Order<double>>& expected_kpoints,
    const int expected_rows,
    const int expected_columns,
    const double tolerance) const
{
    if (kind_ != expected_kind || spin_count <= 0 ||
        expected_kpoints.empty() || expected_rows <= 0 ||
        expected_columns <= 0 || !(tolerance > 0.0) ||
        !std::isfinite(tolerance) ||
        entries_.size() !=
            static_cast<std::size_t>(spin_count) * expected_kpoints.size())
    {
        throw std::invalid_argument(
            "QSGW Vxc manifest does not match the requested dataset");
    }
    for (int spin = 0; spin < spin_count; ++spin)
    {
        for (std::size_t kpoint = 0; kpoint < expected_kpoints.size(); ++kpoint)
        {
            const VxcManifestEntry& entry =
                at(spin, static_cast<int>(kpoint));
            const Vector3_Order<double>& expected = expected_kpoints[kpoint];
            if (entry.rows != expected_rows ||
                entry.columns != expected_columns ||
                periodic_distance(entry.kpoint.x, expected.x) > tolerance ||
                periodic_distance(entry.kpoint.y, expected.y) > tolerance ||
                periodic_distance(entry.kpoint.z, expected.z) > tolerance)
            {
                throw std::invalid_argument(
                    "QSGW Vxc manifest k coordinate does not match its dataset index");
            }
        }
    }
}

void VxcManifest::validate_file_hashes(
    const std::string& base_directory) const
{
    if (base_directory.empty())
    {
        throw std::invalid_argument(
            "QSGW Vxc manifest base directory must not be empty");
    }
    const std::filesystem::path base(base_directory);
    for (const auto& item : entries_)
    {
        const VxcManifestEntry& entry = item.second;
        const std::filesystem::path path = base / entry.file;
        if (sha256_file(path.string()) != entry.sha256)
        {
            throw std::invalid_argument(
                "QSGW Vxc input SHA256 mismatch: " + path.string());
        }
    }
}

const VxcManifestEntry& VxcManifest::at(const int spin,
                                        const int k_index) const
{
    const auto entry = entries_.find({spin, k_index});
    if (entry == entries_.end())
    {
        throw std::out_of_range(
            "QSGW Vxc manifest does not contain the requested spin/k entry");
    }
    return entry->second;
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
            std::string hash;
            std::string key;
            if (!(fields >> hash >> key >> rows))
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
            std::string hash;
            std::string key;
            if (!(fields >> hash >> key >> columns))
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

Matz project_vxc_nao_to_fixed_basis(const Matz& vxc_nao,
                                    const MeanField& reference,
                                    const int spin,
                                    const int kpoint)
{
    validate_reference_wfc(reference, spin, kpoint);
    validate_hermitian_matrix(vxc_nao, reference.get_n_aos(),
                              "QSGW NAO Vxc matrix");
    Matz result(reference.get_n_bands(), reference.get_n_bands(),
                MAJOR::ROW);
    for (int bra = 0; bra < reference.get_n_bands(); ++bra)
    {
        for (int ket = 0; ket < reference.get_n_bands(); ++ket)
        {
            for (int spinor = 0; spinor < reference.get_n_spinor(); ++spinor)
            {
                const ComplexMatrix& wfc =
                    reference.get_eigenvectors()
                        .at(spin)
                        .at(spinor)
                        .at(kpoint);
                for (int row = 0; row < reference.get_n_aos(); ++row)
                {
                    for (int column = 0; column < reference.get_n_aos();
                         ++column)
                    {
                        result(bra, ket) +=
                            std::conj(wfc(bra, row)) * vxc_nao(row, column) *
                            wfc(ket, column);
                    }
                }
            }
        }
    }
    validate_hermitian_matrix(result, reference.get_n_bands(),
                              "QSGW projected Vxc matrix");
    return result;
}

Matz prepare_vxc_in_fixed_state_basis(const Matz& input,
                                      const VxcBasis basis,
                                      const MeanField& reference,
                                      const int spin,
                                      const int kpoint)
{
    validate_reference_wfc(reference, spin, kpoint);
    if (basis == VxcBasis::Nao)
    {
        return project_vxc_nao_to_fixed_basis(
            input, reference, spin, kpoint);
    }
    if (basis == VxcBasis::State)
    {
        validate_hermitian_matrix(input, reference.get_n_bands(),
                                  "QSGW state-basis Vxc matrix");
        return input.copy();
    }
    throw std::invalid_argument("QSGW Vxc basis is invalid");
}

} // namespace qsgw
} // namespace librpa_int
