#pragma once

#include "../core/orbital_weights.h"
#include "librpa_orbital.h"

namespace librpa_int::api
{
LibrpaOrbitalResult* own_orbital_result(OrbitalWeightResult values);
//! The wrapper owns no numerical arrays; the supplied result must outlive it.
LibrpaOrbitalResult* borrow_orbital_result(const OrbitalWeightResult& values);
}  // namespace librpa_int::api
