/*!
 * Common (small) helper functions to initialize dataset components from parsed options.
 */
// Public API headers
#include "librpa_enums.h"
#include "librpa_options.h"

// Internal headers
#include "dataset.h"

namespace librpa_int
{

void initialize_ds_tfgrids(Dataset &ds, const LibrpaOptions &opts);

void initialize_ds_atpairs_local(Dataset &ds, LibrpaParallelRouting routing);

// Initialize response function component
void initialize_ds_chi0(Dataset &ds, const LibrpaOptions &opts) noexcept;

// Initialize EXX component
void initialize_ds_exx(Dataset &ds, const LibrpaOptions &opts) noexcept;

// Initialize G0W0 component
void initialize_ds_g0w0(Dataset &ds, const LibrpaOptions &opts) noexcept;

}
