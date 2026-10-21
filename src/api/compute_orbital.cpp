#include <algorithm>
#include <memory>
#include <string>

#include "../core/orbital_weights.h"
#include "../utils/error.h"
#include "instance_manager.h"
#include "librpa_orbital.h"

struct LibrpaOrbitalResult
{
    librpa_int::OrbitalWeightResult values;
};

LibrpaOrbitalResult* librpa_analyze_orbitals(const LibrpaHandler* h, int nk, int nao,
                                             const double* packed, double tolerance)
{
    using namespace librpa_int;
    const auto ds = api::get_dataset_instance(h);
    auto result = std::make_unique<LibrpaOrbitalResult>();
    std::string error;
    try
    {
        if (!packed || nk != ds->mf.get_n_kpoints() || nao != ds->mf.get_n_aos())
            throw LIBRPA_INVALID_ARGUMENT(
                "orbital-analysis overlap dimensions do not match KS data");
        std::vector<ComplexMatrix> overlap;
        const auto block = static_cast<std::size_t>(nao) * nao;
        for (int k = 0; k < nk; ++k)
        {
            overlap.emplace_back(nao, nao);
            for (std::size_t i = 0; i < block; ++i)
            {
                const auto offset = 2 * (k * block + i);
                overlap.back().c[i] = {packed[offset], packed[offset + 1]};
            }
        }
        std::vector<int> types(ds->basis_wfc.n_atoms);
        for (std::size_t atom = 0; atom < types.size(); ++atom)
            types[atom] = ds->atoms.types.at(atom);
        result->values = analyze_orbital_weights(ds->basis_wfc, ds->basis_convention, types, ds->mf,
                                                 overlap, ds->pbc.weight_k, tolerance);
    }
    catch (const std::exception& failure)
    {
        error = failure.what();
    }
    const int failed = !error.empty();
    int any_failed = 0;
    ds->comm_h.allreduce(&failed, &any_failed, 1, MPI_MAX);
    if (any_failed)
        throw LIBRPA_RUNTIME_ERROR("orbital analysis: " +
                                   (error.empty() ? "another MPI rank failed" : error));
    return result.release();
}

void librpa_delete_orbital_result(LibrpaOrbitalResult* result) { delete result; }

void librpa_orbital_info(const LibrpaOrbitalResult* result, LibrpaOrbitalInfo* info)
{
    if (!result || !info) throw LIBRPA_INVALID_ARGUMENT("null orbital result/info");
    const auto& r = result->values;
    *info = {r.n_spins,
             r.n_spinor,
             r.n_kpoints,
             r.n_bands,
             r.n_aos,
             static_cast<int>(r.groups.size()),
             r.fermi,
             r.highest_occupied,
             r.min_overlap_eigenvalue,
             r.max_norm_error,
             r.max_sum_error,
             r.max_sqrt_residual};
}

void librpa_orbital_group(const LibrpaOrbitalResult* result, int group, LibrpaOrbitalGroup* view)
{
    if (!result || !view || group < 0) throw LIBRPA_INVALID_ARGUMENT("invalid orbital group query");
    const auto& g = result->values.groups.at(group);
    *view = {g.atom, g.type, g.l, g.channel.c_str()};
}

void librpa_orbital_state(const LibrpaOrbitalResult* result, int spin, int k, int band,
                          LibrpaOrbitalState* view)
{
    if (!result || !view) throw LIBRPA_INVALID_ARGUMENT("null orbital result/state");
    const auto& r = result->values;
    if (spin < 0 || spin >= r.n_spins || k < 0 || k >= r.n_kpoints || band < 0 || band >= r.n_bands)
        throw LIBRPA_INVALID_ARGUMENT("orbital state index out of range");
    const auto& s = r.states[(static_cast<std::size_t>(spin) * r.n_kpoints + k) * r.n_bands + band];
    *view = {s.energy,           s.occupation,      s.kweight,      s.norm, s.mulliken_ao.data(),
             s.lowdin_ao.data(), s.mulliken.data(), s.lowdin.data()};
}

void librpa_orbital_pdos(const LibrpaOrbitalResult* result, int count, const double* grid,
                         double sigma, double* mulliken, double* lowdin)
{
    if (!result || count <= 0 || !grid || !mulliken || !lowdin)
        throw LIBRPA_INVALID_ARGUMENT("invalid orbital PDOS buffers");
    std::vector<double> m, l;
    librpa_int::gaussian_orbital_pdos(result->values, {grid, grid + count}, sigma, m, l);
    std::copy(m.begin(), m.end(), mulliken);
    std::copy(l.begin(), l.end(), lowdin);
}
