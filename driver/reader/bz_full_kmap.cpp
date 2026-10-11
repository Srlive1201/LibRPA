#include "bz_full_kmap.h"

#include <array>
#include <cmath>
#include <istream>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>

#include "../../src/api/instance_manager.h"
#include "../../src/utils/constants.h"
#include "reader_context.h"

namespace librpa::reader
{
bool read_full_kmap(std::istream &input, ReaderContext &ctx, const int nk[3],
                    const std::vector<double> &kvecs, const std::vector<double> &weights,
                    const std::vector<int> &qmap)
{
    using namespace librpa_int;
    std::string tag;
    if (!(input >> tag)) return false;  // Existing files have no optional tail.
    // Preserve the old reader's treatment of unrelated trailing data, including
    // numeric IBZ weight summaries. Only the explicit marker selects this path.
    if (tag != "full_kmap") return false;
    long long count;
    const long long expected = static_cast<long long>(nk[0]) * nk[1] * nk[2];
    if (!(input >> count) || count <= 0 || count != expected ||
        count > std::numeric_limits<int>::max())
        throw std::runtime_error("Invalid full_kmap header or full-grid count");
    const int n_scf = static_cast<int>(weights.size());
    if (n_scf < 1 || qmap.size() != weights.size() || kvecs.size() != 3 * weights.size())
        throw std::runtime_error("Invalid full_kmap SCF buffers");
    for (int ik = 0; ik < n_scf; ++ik)
        if (qmap[ik] != ik)
            throw std::runtime_error("full_kmap requires one Coulomb representative per SCF point");

    auto ds = api::get_dataset_instance(ctx.h);
    auto &pbc = ds->pbc;
    std::vector<std::vector<Vector3_Order<double>>> stars(n_scf);
    std::vector<double> full_kvecs;
    std::vector<int> full_to_scf;
    full_kvecs.reserve(3 * static_cast<std::size_t>(count));
    full_to_scf.reserve(static_cast<std::size_t>(count));
    std::vector<std::set<std::array<int, 3>>> star_keys(n_scf);
    std::set<std::array<int, 3>> seen;
    std::array<double, 3> origin;
    // Validate a uniform mesh modulo reciprocal lattice vectors, including its
    // actual offset. Do not infer the producer's IBZ from a regenerated mesh.
    const auto grid_key = [&](const Vector3_Order<double> &k)
    {
        const auto frac = pbc.latvec * k;
        const double f[3] = {frac.x, frac.y, frac.z};
        std::array<int, 3> key;
        for (int axis = 0; axis < 3; ++axis)
        {
            const double x = (f[axis] - origin[axis]) * nk[axis];
            if (!std::isfinite(x) || std::abs(x - std::round(x)) > 1e-5)
                throw std::runtime_error("full_kmap point is outside the declared uniform grid");
            const double reduced = std::fmod(std::round(x), nk[axis]);
            key[axis] = (static_cast<int>(reduced) + nk[axis]) % nk[axis];
        }
        return key;
    };
    for (int row = 0; row < count; ++row)
    {
        int id, rep;
        Vector3_Order<double> k;
        if (!(input >> id >> k.x >> k.y >> k.z >> rep) || id != row + 1 || rep < 1 || rep > n_scf ||
            !std::isfinite(k.x) || !std::isfinite(k.y) || !std::isfinite(k.z))
            throw std::runtime_error("Invalid full_kmap row or SCF representative index");
        full_kvecs.insert(full_kvecs.end(), {k.x, k.y, k.z});
        full_to_scf.push_back(rep - 1);
        k /= TWO_PI;  // Same internal units as PeriodicBoundaryData::klist.
        if (row == 0)
        {
            const auto frac = pbc.latvec * k;
            origin = {frac.x, frac.y, frac.z};
        }
        const auto key = grid_key(k);
        if (!seen.insert(key).second)
            throw std::runtime_error("Duplicate full_kmap point modulo reciprocal lattice vectors");
        stars[rep - 1].emplace_back(k);
        star_keys[rep - 1].insert(key);
    }
    if (input >> tag) throw std::runtime_error("Unexpected data after full_kmap");
    for (int ik = 0; ik < n_scf; ++ik)
    {
        const Vector3_Order<double> k{kvecs[3 * ik] / TWO_PI, kvecs[3 * ik + 1] / TWO_PI,
                                      kvecs[3 * ik + 2] / TWO_PI};
        if (star_keys[ik].count(grid_key(k)) != 1)
            throw std::runtime_error("SCF representative is absent from its full_kmap star");
        const double weight = static_cast<double>(stars[ik].size()) / count;
        if (std::abs(weight - weights[ik]) > 1e-6)
            throw std::runtime_error("Folded SCF weight disagrees with full_kmap multiplicity");
    }
    ctx.h.set_irreducible_kgrids_kvec(nk[0], nk[1], nk[2], kvecs, full_kvecs, full_to_scf);
    return true;
}
}  // namespace librpa::reader
