#pragma once
#include "../core/meanfield.h"
#include "../math/matrix_m.h"

namespace librpa_int
{
namespace qsgw
{
enum class VxcBasis
{
    Nao,
    State,
};

Matz project_vxc_nao_to_fixed_basis(const Matz& vxc_nao,
                                    const MeanField& reference,
                                    int spin,
                                    int kpoint);

Matz prepare_vxc_in_fixed_state_basis(const Matz& input,
                                      VxcBasis basis,
                                      const MeanField& reference,
                                      int spin,
                                      int kpoint);

} // namespace qsgw
} // namespace librpa_int
