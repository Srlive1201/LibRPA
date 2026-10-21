#pragma once

#include "../../src/core/qsgw/matrix_map.h"
#include "../../src/core/qsgw/vxc_projection.h"

#include <string>

namespace librpa_int
{
namespace qsgw
{

// Read full Vxc matrices in the same spin/k order as the SCF or band reference.
// An empty prefix selects the producer's standard SCF/band filename prefix.
// ABACUS out_mat_xc is in the full KS-state basis, including for SOC, and
// uses the s4 filename suffix for its single two-component spinor channel.
SpinKMatrixMap read_qsgw_vxc(const std::string& input_directory,
                            const std::string& prefix,
                            const MeanField& reference,
                            bool aims_input,
                            bool band_path,
                            VxcBasis basis);

} // namespace qsgw
} // namespace librpa_int
