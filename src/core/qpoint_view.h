/*!
 * @file qpoint_view.h
 * @brief Derived q-point views for response and screened-Coulomb calculations.
 */
#pragma once

#include <map>
#include <vector>

#include "pbc.h"
#include "qpoint_types.h"
#include "symmetry_context.h"

namespace librpa_int
{

SymmetryQPointView build_symmetry_qpoint_view(
    const SymmetryContext& ctx,
    const PeriodicBoundaryData& pbc,
    bool use_symmetry);

} // namespace librpa_int
