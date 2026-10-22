#include "band_remap.h"

#include <string>

#include "../../io/global_io.h"
#include "../../io/stl_io_helper.h"
#include "../../utils/error.h"

namespace librpa_int::qsgw
{
AtomPairBvKRemap<atom_t> build_qsgw_band_bvk_remap(const Atoms &atoms, const PeriodicBoundaryData &pbc,
                                              const int remap_convention)
{
    using global::ofs_myid;

    if (remap_convention != 0 && remap_convention != 1)
        throw LIBRPA_RUNTIME_ERROR("Invalid BvK remap convention: " +
                                   std::to_string(remap_convention));

    const bool is_ready =
        pbc.is_latt_set() && atoms.is_set() && atoms.is_frac_set() && !pbc.Rlist.empty();
    AtomPairBvKRemap<atom_t> remap;
    if (is_ready)
    {
        remap.build(atoms.coords_frac, pbc.Rlist, pbc.period, pbc.latvec, remap_convention);
    }

    ofs_myid << "Final BvK remapping for band interpolation:\n";
    ofs_myid << "| remap convention  : " << remap_convention << "\n";
    ofs_myid << "| input ready       : " << (is_ready ? "true" : "false") << "\n";
    ofs_myid << "| lattice set       : " << (pbc.is_latt_set() ? "true" : "false") << "\n";
    ofs_myid << "| atoms set         : " << (atoms.is_set() ? "true" : "false") << "\n";
    ofs_myid << "| atom frac set     : " << (atoms.is_frac_set() ? "true" : "false") << "\n";
    ofs_myid << "| period            : " << pbc.period << "\n";
    ofs_myid << "| Rlist size        : " << pbc.Rlist.size() << "\n";
    ofs_myid << "| atom-pair entries : " << remap.size() << "\n";
    for (const auto &[atom_pair, R_map] : remap.data())
    {
        for (const auto &[R, R_bvks] : R_map)
        {
            // Only print 1-to-1 remapping with debug
            if (R_bvks.size() < 2 && !global::should_output(LIBRPA_VERBOSE_DEBUG)) continue;
            ofs_myid << "| " << atom_pair << " " << R << " -> " << R_bvks << "\n";
        }
    }
    ofs_myid << std::flush;
    return remap;
}
}  // namespace librpa_int::qsgw
