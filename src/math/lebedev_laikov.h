#pragma once

#include <vector>

namespace librpa_int::lebedev_laikov {

struct Grid
{
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;
    std::vector<double> weights;
};

std::vector<int> available_orders();

bool is_available_order(int order);

Grid grid(int order);

} // namespace librpa_int::lebedev_laikov
