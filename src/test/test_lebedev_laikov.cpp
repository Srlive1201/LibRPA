#include "../math/lebedev_laikov.h"

#include <cassert>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

using namespace librpa_int;

static void check_grid(const lebedev_laikov::Grid &grid, std::size_t expected_order)
{
    assert(grid.x.size() == expected_order);
    assert(grid.y.size() == expected_order);
    assert(grid.z.size() == expected_order);
    assert(grid.weights.size() == expected_order);

    double sum_w = 0.0;
    double sum_x = 0.0;
    double sum_y = 0.0;
    double sum_z = 0.0;
    double sum_x2 = 0.0;
    double sum_y2 = 0.0;
    double sum_z2 = 0.0;
    double sum_xy = 0.0;

    for (std::size_t i = 0; i != grid.weights.size(); ++i)
    {
        const double x = grid.x[i];
        const double y = grid.y[i];
        const double z = grid.z[i];
        const double w = grid.weights[i];
        const double r2 = x * x + y * y + z * z;

        assert(std::abs(r2 - 1.0) < 2.0e-15);

        sum_w += w;
        sum_x += w * x;
        sum_y += w * y;
        sum_z += w * z;
        sum_x2 += w * x * x;
        sum_y2 += w * y * y;
        sum_z2 += w * z * z;
        sum_xy += w * x * y;
    }

    assert(std::abs(sum_w - 1.0) < 5.0e-14);
    assert(std::abs(sum_x) < 2.0e-16);
    assert(std::abs(sum_y) < 2.0e-16);
    assert(std::abs(sum_z) < 2.0e-16);
    assert(std::abs(sum_x2 - 1.0 / 3.0) < 2.0e-14);
    assert(std::abs(sum_y2 - 1.0 / 3.0) < 2.0e-14);
    assert(std::abs(sum_z2 - 1.0 / 3.0) < 2.0e-14);
    assert(std::abs(sum_xy) < 2.0e-16);
}

static void check_grid_5810_high_moments()
{
    const auto grid = lebedev_laikov::grid(5810);
    double sum_x4 = 0.0;
    double sum_x2y2 = 0.0;

    for (std::size_t i = 0; i != grid.weights.size(); ++i)
    {
        const double x = grid.x[i];
        const double y = grid.y[i];
        const double w = grid.weights[i];
        sum_x4 += w * x * x * x * x;
        sum_x2y2 += w * x * x * y * y;
    }

    assert(std::abs(sum_x4 - 1.0 / 5.0) < 2.0e-14);
    assert(std::abs(sum_x2y2 - 1.0 / 15.0) < 2.0e-14);
}

static void test_available_orders()
{
    const std::vector<int> expected = {6,    14,   26,   38,   50,   74,   86,   110,  146,  170,  194,
                                       230,  266,  302,  350,  434,  590,  770,  974,  1202, 1454, 1730,
                                       2030, 2354, 2702, 3074, 3470, 3890, 4334, 4802, 5294, 5810};
    assert(lebedev_laikov::available_orders() == expected);

    for (const int order : expected)
    {
        assert(lebedev_laikov::is_available_order(order));
        check_grid(lebedev_laikov::grid(order), static_cast<std::size_t>(order));
    }

    assert(!lebedev_laikov::is_available_order(42));
    try
    {
        (void)lebedev_laikov::grid(42);
        assert(false);
    }
    catch (const std::runtime_error &)
    {
    }
}

int main()
{
    test_available_orders();
    check_grid_5810_high_moments();
    return 0;
}
