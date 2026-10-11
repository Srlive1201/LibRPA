#pragma once

#include <iosfwd>
#include <vector>

#include "reader_context.h"

namespace librpa::reader
{
bool read_full_kmap(std::istream &input, ReaderContext &ctx, const int nk[3],
                    const std::vector<double> &kvecs, const std::vector<double> &weights,
                    const std::vector<int> &qmap);
}  // namespace librpa::reader
