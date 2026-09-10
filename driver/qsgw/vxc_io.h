#pragma once

#include "../../src/qsgw/vxc_projection.h"
#include "../../src/math/matrix_m.h"

#include <istream>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace librpa_int
{
namespace qsgw
{

enum class VxcDatasetKind
{
    ScfGrid,
    BandPath,
};

enum class VxcUnits
{
    Hartree,
    Rydberg,
};

enum class VxcGauge
{
    AoBloch,
    Mf0State,
};

struct VxcManifestEntry
{
    int spin = -1;
    int k_index = -1;
    Vector3_Order<double> kpoint;
    int rows = -1;
    int columns = -1;
    std::string file;
};

class VxcManifest
{
public:
    static VxcManifest parse(std::istream& input,
                             const std::string& source_name);

    void validate(VxcDatasetKind expected_kind,
                  int spin_count,
                  const std::vector<Vector3_Order<double>>& expected_kpoints,
                  int expected_rows,
                  int expected_columns,
                  double tolerance) const;


    const std::string& producer() const noexcept { return producer_; }
    VxcUnits units() const noexcept { return units_; }
    VxcBasis basis() const noexcept { return basis_; }
    VxcGauge gauge() const noexcept { return gauge_; }
    VxcDatasetKind kind() const noexcept { return kind_; }
    const VxcManifestEntry& at(int spin, int k_index) const;

private:
    VxcDatasetKind kind_ = VxcDatasetKind::ScfGrid;
    std::string producer_;
    VxcUnits units_ = VxcUnits::Hartree;
    VxcBasis basis_ = VxcBasis::State;
    VxcGauge gauge_ = VxcGauge::Mf0State;
    std::map<std::pair<int, int>, VxcManifestEntry> entries_;
};

Matz read_abacus_vxc_ha(std::istream& input,
                        const std::string& source_name);

} // namespace qsgw
} // namespace librpa_int
