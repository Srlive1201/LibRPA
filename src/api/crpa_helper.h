#pragma once

#include "../core/crpa_workflow.h"
#include "dataset.h"
#include "librpa_crpa.h"

namespace librpa_int
{
namespace api
{
enum class CrpaApiMode
{
    Standalone,
    SharedGw,
    ResponseOnly
};
CrpaPreparedModel prepare_crpa_api_model(Dataset& dataset, const LibrpaOptions& options,
                                         const LibrpaCrpaInput* input,
                                         CrpaApiMode mode = CrpaApiMode::Standalone);

//! Consume the prepared model while preserving the caller's resident response.
LibrpaCrpaResult* evaluate_crpa_api_from_response(Dataset& dataset, const LibrpaOptions& options,
                                                  CrpaPreparedModel& model,
                                                  const Chi0QMap& full_response);
}  // namespace api
}  // namespace librpa_int
