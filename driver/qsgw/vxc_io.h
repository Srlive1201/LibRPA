#pragma once

#include "../../src/qsgw/matrix_map.h"
#include "../../src/qsgw/vxc_projection.h"

#include <string>

namespace librpa_int
{
namespace qsgw
{

// Read full Vxc matrices in the same spin/k order as the SCF or band reference.
// An empty prefix selects the producer's standard SCF/band filename prefix.
SpinKMatrixMap read_qsgw_vxc(const std::string& input_directory,
                            const std::string& prefix,
                            const MeanField& reference,
                            bool aims_input,
                            bool band_path,
                            VxcBasis basis);

} // namespace qsgw
} // namespace librpa_int
