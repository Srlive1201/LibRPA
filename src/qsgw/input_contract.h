#pragma once

#include <filesystem>
#include <istream>
#include <map>
#include <string>
#include <vector>

namespace librpa_int
{
namespace qsgw
{

enum class QsgwProducer
{
    Abacus,
    FhiAims,
};

enum class HeadwingGridMode
{
    Disabled,
    ScfGrid,
};

enum class HeadwingUpdateMode
{
    None,
    FixedReference,
};

enum class BandUpdateMode
{
    Off,
    FixedBasisRotation,
};

struct QsgwInputFile
{
    std::string role;
    std::string sha256;
    std::string file;
};

struct BandReferencePaths
{
    std::filesystem::path kpoints;
    std::vector<std::filesystem::path> eigenvalues;
    std::vector<std::filesystem::path> wavefunctions;
};

BandReferencePaths resolve_band_reference_paths(
    const std::string& input_directory,
    const std::string& band_kpath_filename,
    int n_kpoints);

std::vector<std::filesystem::path> resolve_same_grid_velocity_paths(
    QsgwProducer producer,
    const std::string& input_directory,
    int n_kpoints);

class QsgwInputContract
{
public:
    static QsgwInputContract parse(std::istream& input,
                                   const std::string& source_name);

    void validate_file_hashes(const std::string& base_directory) const;
    const std::vector<QsgwInputFile>& files(const std::string& role) const;

    QsgwProducer producer() const noexcept { return producer_; }
    HeadwingGridMode headwing_grid() const noexcept { return headwing_grid_; }
    HeadwingUpdateMode headwing_update() const noexcept
    {
        return headwing_update_;
    }
    BandUpdateMode band_update() const noexcept { return band_update_; }
    int n_spins() const noexcept { return n_spins_; }
    int n_bands() const noexcept { return n_bands_; }
    int n_aos() const noexcept { return n_aos_; }
    int n_scf_kpoints() const noexcept { return n_scf_kpoints_; }
    int n_headwing_kpoints() const noexcept { return n_headwing_kpoints_; }
    int n_band_kpoints() const noexcept { return n_band_kpoints_; }

private:
    QsgwProducer producer_ = QsgwProducer::Abacus;
    HeadwingGridMode headwing_grid_ = HeadwingGridMode::Disabled;
    HeadwingUpdateMode headwing_update_ = HeadwingUpdateMode::None;
    BandUpdateMode band_update_ = BandUpdateMode::Off;
    int n_spins_ = 0;
    int n_bands_ = 0;
    int n_aos_ = 0;
    int n_scf_kpoints_ = 0;
    int n_headwing_kpoints_ = 0;
    int n_band_kpoints_ = 0;
    std::map<std::string, std::vector<QsgwInputFile>> files_;
};

void validate_scf_input_binding(
    const QsgwInputContract& contract,
    const std::string& contract_base_directory,
    const std::filesystem::path& eigenvalue_file,
    const std::vector<std::filesystem::path>& wavefunction_files,
    const std::filesystem::path& kpoint_file,
    const std::vector<std::filesystem::path>& reader_static_files);

void validate_band_reference_binding(
    const QsgwInputContract& contract,
    const std::string& contract_base_directory,
    const std::string& input_directory,
    const std::string& band_kpath_filename);

void validate_qsgw_execution_modes(const QsgwInputContract& contract,
                                   HeadwingGridMode headwing_grid,
                                   bool update_band);

} // namespace qsgw
} // namespace librpa_int
