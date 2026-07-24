#include "input_contract.h"
#include "sha256.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <set>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace librpa_int
{
namespace qsgw
{
namespace
{

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

int parse_integer(const std::string& value,
                  const std::string& key,
                  const std::string& source_name)
{
    try
    {
        std::size_t consumed = 0;
        const int result = std::stoi(value, &consumed);
        if (consumed != value.size()) throw std::invalid_argument("suffix");
        return result;
    }
    catch (const std::exception&)
    {
        throw std::invalid_argument(
            "Invalid QSGW input-contract integer " + key + " in " +
            source_name);
    }
}

int read_pyatb_kpoint_count(const std::filesystem::path& kpoint_file,
                            const int active_kpoints)
{
    std::ifstream input(kpoint_file);
    int n_basis = 0;
    int n_states = 0;
    int n_spins = 0;
    int source_kpoints = 0;
    if (!(input >> n_basis >> n_states >> n_spins >> source_kpoints) ||
        n_basis <= 0 || n_states <= 0 || n_spins <= 0 ||
        source_kpoints < active_kpoints)
    {
        throw std::invalid_argument(
            "QSGW cannot determine complete PyATB head-only k-point coverage from " +
            kpoint_file.string());
    }
    return source_kpoints;
}

bool safe_relative_path(const std::string& value)
{
    if (value.empty() || value.find('\\') != std::string::npos) return false;
    const std::filesystem::path path(value);
    if (path.is_absolute()) return false;
    for (const auto& component : path)
    {
        if (component == "..") return false;
    }
    return true;
}

const std::set<std::string>& allowed_roles()
{
    static const std::set<std::string> roles{
        "mf0_eigenvalues",
        "mf0_wavefunctions",
        "scf_kpoints",
        "vxc_scf_manifest",
        "reader_static",
        "velocity_mf0",
        "band_mf0_eigenvalues",
        "band_mf0_wavefunctions",
        "band_kpoints",
        "vxc_band_manifest",
    };
    return roles;
}

bool has_role(const QsgwInputContract& contract, const std::string& role)
{
    return !contract.files(role).empty();
}

void require_roles(const QsgwInputContract& contract,
                   const std::initializer_list<const char*> roles,
                   const std::string& source_name)
{
    for (const char* role : roles)
    {
        if (!has_role(contract, role))
        {
            throw std::invalid_argument(
                "Missing QSGW input-contract role " + std::string(role) +
                " in " + source_name);
        }
    }
}

void reject_roles(const QsgwInputContract& contract,
                  const std::initializer_list<const char*> roles,
                  const std::string& source_name)
{
    for (const char* role : roles)
    {
        if (has_role(contract, role))
        {
            throw std::invalid_argument(
                "Unexpected QSGW input-contract role " + std::string(role) +
                " in " + source_name);
        }
    }
}

} // namespace

BandReferencePaths resolve_band_reference_paths(
    const std::string& input_directory,
    const std::string& band_kpath_filename,
    const int n_kpoints)
{
    if (input_directory.empty() || band_kpath_filename.empty() ||
        n_kpoints <= 0)
    {
        throw std::invalid_argument(
            "QSGW band reference input needs a directory, k-path file, and positive k-point count");
    }

    const std::filesystem::path input(input_directory);
    BandReferencePaths paths;
    paths.kpoints = std::filesystem::absolute(
                        input / band_kpath_filename)
                        .lexically_normal();
    paths.eigenvalues.reserve(static_cast<std::size_t>(n_kpoints));
    paths.wavefunctions.reserve(static_cast<std::size_t>(n_kpoints));
    for (int kpoint = 1; kpoint <= n_kpoints; ++kpoint)
    {
        std::ostringstream index;
        index << std::setw(5) << std::setfill('0') << kpoint;
        paths.eigenvalues.push_back(
            std::filesystem::absolute(
                input /
                ("band_KS_eigenvalue_k_" + index.str() + ".txt"))
                .lexically_normal());
        paths.wavefunctions.push_back(
            std::filesystem::absolute(
                input /
                ("band_KS_eigenvector_k_" + index.str() + ".txt"))
                .lexically_normal());
    }
    return paths;
}

std::vector<std::filesystem::path> resolve_same_grid_velocity_paths(
    const QsgwProducer producer,
    const std::string& input_directory,
    const int n_kpoints)
{
    if (input_directory.empty() || n_kpoints <= 0)
    {
        throw std::invalid_argument(
            "QSGW same-grid velocity input needs a directory and positive k-point count");
    }

    const std::filesystem::path input(input_directory);
    std::vector<std::filesystem::path> paths;
    if (producer == QsgwProducer::FhiAims)
    {
        paths.reserve(static_cast<std::size_t>(n_kpoints));
        for (int kpoint = 1; kpoint <= n_kpoints; ++kpoint)
        {
            std::ostringstream filename;
            filename << "mommat_ks_kpt_" << std::setw(6)
                     << std::setfill('0') << kpoint << ".dat";
            paths.push_back(
                std::filesystem::absolute(input / filename.str())
                    .lexically_normal());
        }
    }
    else
    {
        const std::filesystem::path pyatb_velocity =
            input / "pyatb_librpa_df" / "velocity_matrix";
        if (std::filesystem::exists(pyatb_velocity))
        {
            const std::filesystem::path pyatb =
                input / "pyatb_librpa_df";
            const std::filesystem::path pyatb_kpoints =
                pyatb / "k_path_info";
            const int source_kpoints =
                read_pyatb_kpoint_count(pyatb_kpoints, n_kpoints);
            paths.reserve(static_cast<std::size_t>(source_kpoints) + 3);
            paths.push_back(pyatb_kpoints);
            paths.push_back(pyatb / "band_out");
            for (int kpoint = 0; kpoint < source_kpoints; ++kpoint)
            {
                paths.push_back(
                    pyatb /
                    ("KS_eigenvector_" + std::to_string(kpoint) + ".dat"));
            }
            paths.push_back(pyatb_velocity);
            for (std::filesystem::path& path : paths)
            {
                path = std::filesystem::absolute(path).lexically_normal();
            }
        }
        else
        {
            paths.push_back(
                std::filesystem::absolute(input / "velocity_matrix")
                    .lexically_normal());
        }
    }
    return paths;
}

QsgwInputContract QsgwInputContract::parse(
    std::istream& input,
    const std::string& source_name)
{
    QsgwInputContract result;
    std::map<std::string, std::string> metadata;
    std::set<std::pair<std::string, std::string>> role_files;
    bool saw_magic = false;
    bool saw_table_header = false;
    std::string line;
    while (std::getline(input, line))
    {
        const std::string content = trim(line);
        if (content.empty()) continue;
        if (!saw_magic)
        {
            if (content != "# librpa-qsgw-input-contract-v1")
            {
                throw std::invalid_argument(
                    "Invalid QSGW input-contract header in " + source_name);
            }
            saw_magic = true;
            continue;
        }
        if (!saw_table_header)
        {
            std::istringstream fields(content);
            std::string key;
            std::string value;
            fields >> key;
            if (lowercase(key) == "role")
            {
                std::string sha256;
                std::string file;
                std::string extra;
                if (!(fields >> sha256 >> file) || fields >> extra ||
                    lowercase(sha256) != "sha256" ||
                    lowercase(file) != "file")
                {
                    throw std::invalid_argument(
                        "Invalid QSGW input-contract table header in " +
                        source_name);
                }
                saw_table_header = true;
                continue;
            }
            std::string extra;
            if (key.empty() || !(fields >> value) || fields >> extra)
            {
                throw std::invalid_argument(
                    "Malformed QSGW input-contract metadata in " +
                    source_name);
            }
            key = lowercase(key);
            if (!metadata.emplace(key, value).second)
            {
                throw std::invalid_argument(
                    "Duplicate QSGW input-contract metadata in " +
                    source_name);
            }
            continue;
        }

        QsgwInputFile file;
        std::istringstream fields(content);
        std::string extra;
        if (!(fields >> file.role >> file.sha256 >> file.file) ||
            fields >> extra)
        {
            throw std::invalid_argument(
                "Malformed QSGW input-contract file entry in " +
                source_name);
        }
        file.role = lowercase(file.role);
        if (allowed_roles().count(file.role) == 0 ||
            !is_sha256_hex(file.sha256) || !safe_relative_path(file.file) ||
            !role_files.emplace(file.role, file.file).second)
        {
            throw std::invalid_argument(
                "Invalid QSGW input-contract file entry in " + source_name);
        }
        result.files_[file.role].push_back(std::move(file));
    }

    const std::set<std::string> required_metadata{
        "producer", "internal_energy_units", "mf0_basis", "mf0_gauge",
        "n_spins", "n_bands", "n_aos", "n_scf_kpoints",
        "n_headwing_kpoints", "n_band_kpoints", "headwing_grid",
        "headwing_update", "hartree_update", "band_update",
    };
    if (!saw_magic || !saw_table_header || result.files_.empty() ||
        metadata.size() != required_metadata.size())
    {
        throw std::invalid_argument(
            "Incomplete QSGW input contract " + source_name);
    }
    for (const std::string& key : required_metadata)
    {
        if (metadata.count(key) == 0)
        {
            throw std::invalid_argument(
                "Missing QSGW input-contract metadata " + key + " in " +
                source_name);
        }
    }

    const std::string producer = lowercase(metadata.at("producer"));
    if (producer == "abacus")
        result.producer_ = QsgwProducer::Abacus;
    else if (producer == "fhi-aims")
        result.producer_ = QsgwProducer::FhiAims;
    else
        throw std::invalid_argument(
            "Unsupported QSGW input-contract producer in " + source_name);

    if (lowercase(metadata.at("internal_energy_units")) != "hartree" ||
        lowercase(metadata.at("mf0_basis")) !=
            "state_coefficients_in_nao" ||
        lowercase(metadata.at("mf0_gauge")) != "producer_state")
    {
        throw std::invalid_argument(
            "Invalid QSGW mf0 units, basis, or gauge in " + source_name);
    }

    result.n_spins_ = parse_integer(
        metadata.at("n_spins"), "n_spins", source_name);
    result.n_bands_ = parse_integer(
        metadata.at("n_bands"), "n_bands", source_name);
    result.n_aos_ = parse_integer(
        metadata.at("n_aos"), "n_aos", source_name);
    result.n_scf_kpoints_ = parse_integer(
        metadata.at("n_scf_kpoints"), "n_scf_kpoints", source_name);
    result.n_headwing_kpoints_ = parse_integer(
        metadata.at("n_headwing_kpoints"), "n_headwing_kpoints", source_name);
    result.n_band_kpoints_ = parse_integer(
        metadata.at("n_band_kpoints"), "n_band_kpoints", source_name);
    if (result.n_spins_ <= 0 || result.n_bands_ <= 0 ||
        result.n_aos_ <= 0 || result.n_scf_kpoints_ <= 0 ||
        result.n_headwing_kpoints_ < 0 || result.n_band_kpoints_ < 0)
    {
        throw std::invalid_argument(
            "Invalid QSGW input-contract dimensions in " + source_name);
    }

    const std::string headwing_grid =
        lowercase(metadata.at("headwing_grid"));
    if (headwing_grid == "disabled")
        result.headwing_grid_ = HeadwingGridMode::Disabled;
    else if (headwing_grid == "scf")
        result.headwing_grid_ = HeadwingGridMode::ScfGrid;
    else
        throw std::invalid_argument(
            "QSGW supports only disabled or same-grid head-only input in " +
            source_name);

    const std::string headwing_update =
        lowercase(metadata.at("headwing_update"));
    if (headwing_update == "none")
        result.headwing_update_ = HeadwingUpdateMode::None;
    else if (headwing_update == "fixed_reference")
        result.headwing_update_ = HeadwingUpdateMode::FixedReference;
    else
        throw std::invalid_argument(
            "QSGW supports only fixed-basis same-grid head updates in " +
            source_name);

    const std::string hartree = lowercase(metadata.at("hartree_update"));
    if (hartree != "off")
        throw std::invalid_argument(
            "QSGW Hartree updates are not enabled in " + source_name);

    const std::string band = lowercase(metadata.at("band_update"));
    if (band == "off")
        result.band_update_ = BandUpdateMode::Off;
    else if (band == "fixed_basis_rotation")
        result.band_update_ = BandUpdateMode::FixedBasisRotation;
    else
        throw std::invalid_argument(
            "Unsupported QSGW band update mode in " + source_name);

    require_roles(result,
                  {"mf0_eigenvalues", "mf0_wavefunctions", "scf_kpoints",
                   "vxc_scf_manifest", "reader_static"},
                  source_name);

    if (result.headwing_grid_ == HeadwingGridMode::Disabled)
    {
        if (result.n_headwing_kpoints_ != 0 ||
            result.headwing_update_ != HeadwingUpdateMode::None)
        {
            throw std::invalid_argument(
                "Disabled QSGW head-wing has inconsistent dimensions or update mode in " +
                source_name);
        }
        reject_roles(result, {"velocity_mf0"}, source_name);
    }
    else if (result.headwing_grid_ == HeadwingGridMode::ScfGrid)
    {
        if (result.n_headwing_kpoints_ != result.n_scf_kpoints_ ||
            result.headwing_update_ !=
                HeadwingUpdateMode::FixedReference)
        {
            throw std::invalid_argument(
                "Same-grid QSGW head must retain velocity in the immutable mf0 basis in " +
                source_name);
        }
        require_roles(result, {"velocity_mf0"}, source_name);
    }

    const std::initializer_list<const char*> band_roles{
        "band_mf0_eigenvalues", "band_mf0_wavefunctions", "band_kpoints",
        "vxc_band_manifest"};
    if (result.band_update_ == BandUpdateMode::Off)
    {
        if (result.n_band_kpoints_ != 0)
        {
            throw std::invalid_argument(
                "Disabled QSGW band update has nonzero band k-points in " +
                source_name);
        }
        reject_roles(result, band_roles, source_name);
    }
    else
    {
        if (result.n_band_kpoints_ <= 0)
        {
            throw std::invalid_argument(
                "QSGW fixed-basis band update needs band k-points in " +
                source_name);
        }
        require_roles(result, band_roles, source_name);
    }
    return result;
}

void QsgwInputContract::validate_file_hashes(
    const std::string& base_directory) const
{
    if (base_directory.empty())
    {
        throw std::invalid_argument(
            "QSGW input-contract base directory must not be empty");
    }
    const std::filesystem::path base(base_directory);
    for (const auto& role : files_)
    {
        for (const QsgwInputFile& file : role.second)
        {
            const std::filesystem::path path = base / file.file;
            if (sha256_file(path.string()) != file.sha256)
            {
                throw std::invalid_argument(
                    "QSGW input SHA256 mismatch for role " + file.role +
                    ": " + path.string());
            }
        }
    }
}

const std::vector<QsgwInputFile>& QsgwInputContract::files(
    const std::string& role) const
{
    static const std::vector<QsgwInputFile> empty;
    const auto found = files_.find(lowercase(role));
    return found == files_.end() ? empty : found->second;
}

void validate_scf_input_binding(
    const QsgwInputContract& contract,
    const std::string& contract_base_directory,
    const std::filesystem::path& eigenvalue_file,
    const std::vector<std::filesystem::path>& wavefunction_files,
    const std::filesystem::path& kpoint_file,
    const std::vector<std::filesystem::path>& reader_static_files)
{
    if (contract_base_directory.empty() || eigenvalue_file.empty() ||
        wavefunction_files.empty() || kpoint_file.empty() ||
        reader_static_files.empty())
    {
        throw std::invalid_argument(
            "QSGW SCF binding requires complete reader file sets and a contract base directory");
    }

    const auto declared_paths = [&](const std::string& role) {
        std::vector<std::filesystem::path> paths;
        for (const QsgwInputFile& file : contract.files(role))
        {
            paths.push_back(
                std::filesystem::absolute(
                    std::filesystem::path(contract_base_directory) /
                    file.file)
                    .lexically_normal());
        }
        std::sort(paths.begin(), paths.end());
        return paths;
    };
    const auto require_exact = [&](const std::string& role,
                                   std::vector<std::filesystem::path> paths) {
        for (std::filesystem::path& path : paths)
            path = std::filesystem::absolute(path).lexically_normal();
        std::sort(paths.begin(), paths.end());
        if (declared_paths(role) != paths)
        {
            throw std::invalid_argument(
                "QSGW " + role +
                " contract files do not exactly match the SCF files read by the driver");
        }
    };

    require_exact("mf0_eigenvalues", {eigenvalue_file});
    require_exact("mf0_wavefunctions", wavefunction_files);
    require_exact("scf_kpoints", {kpoint_file});
    require_exact("reader_static", reader_static_files);
}

void validate_band_reference_binding(
    const QsgwInputContract& contract,
    const std::string& contract_base_directory,
    const std::string& input_directory,
    const std::string& band_kpath_filename)
{
    if (contract.band_update() != BandUpdateMode::FixedBasisRotation ||
        contract_base_directory.empty())
    {
        throw std::invalid_argument(
            "QSGW band reference binding requires a fixed-basis rotation contract and contract base directory");
    }
    const BandReferencePaths expected = resolve_band_reference_paths(
        input_directory, band_kpath_filename, contract.n_band_kpoints());

    const auto declared_paths = [&](const std::string& role) {
        std::vector<std::filesystem::path> paths;
        for (const QsgwInputFile& file : contract.files(role))
        {
            paths.push_back(
                std::filesystem::absolute(
                    std::filesystem::path(contract_base_directory) /
                    file.file)
                    .lexically_normal());
        }
        std::sort(paths.begin(), paths.end());
        return paths;
    };
    const auto require_exact = [&](const std::string& role,
                                   std::vector<std::filesystem::path> paths) {
        std::sort(paths.begin(), paths.end());
        if (declared_paths(role) != paths)
        {
            throw std::invalid_argument(
                "QSGW " + role +
                " contract files do not exactly match the band files read by the driver");
        }
    };

    require_exact("band_kpoints", {expected.kpoints});
    require_exact("band_mf0_eigenvalues", expected.eigenvalues);
    require_exact("band_mf0_wavefunctions", expected.wavefunctions);
    if (contract.files("vxc_band_manifest").size() != 1)
    {
        throw std::invalid_argument(
            "QSGW band reference contract requires exactly one Vxc manifest");
    }
}

void validate_qsgw_execution_modes(const QsgwInputContract& contract,
                                   const HeadwingGridMode headwing_grid,
                                   const bool update_band)
{
    HeadwingUpdateMode expected_update = HeadwingUpdateMode::None;
    if (headwing_grid == HeadwingGridMode::ScfGrid)
        expected_update = HeadwingUpdateMode::FixedReference;
    const bool headwing_matches =
        contract.headwing_grid() == headwing_grid &&
        contract.headwing_update() == expected_update;
    if (!headwing_matches)
    {
        throw std::invalid_argument(
            "QSGW execution mode and input contract disagree about the head/wing grid or live-update mode");
    }

    const bool band_matches = update_band
        ? contract.band_update() == BandUpdateMode::FixedBasisRotation
        : contract.band_update() == BandUpdateMode::Off;
    if (!band_matches)
    {
        throw std::invalid_argument(
            "QSGW execution mode and input contract disagree about the band update");
    }
}

} // namespace qsgw
} // namespace librpa_int
