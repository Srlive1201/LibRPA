#ifdef NDEBUG
#undef NDEBUG
#endif

#include "../qsgw/convergence.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

using librpa_int::MeanField;
using librpa_int::qsgw::eigenvalue_snapshot;
using librpa_int::qsgw::max_eigenvalue_change;
using librpa_int::qsgw::qsgw_iteration_converged;

namespace
{

template <typename Function>
void assert_throws(Function&& function)
{
    bool threw = false;
    try
    {
        function();
    }
    catch (const std::exception&)
    {
        threw = true;
    }
    assert(threw);
}

void test_maximum_change_and_minimum_iteration_gate()
{
    MeanField meanfield(2, 2, 3, 1, 1);
    for (int spin = 0; spin < 2; ++spin)
    {
        for (int kpoint = 0; kpoint < 2; ++kpoint)
        {
            for (int band = 0; band < 3; ++band)
            {
                meanfield.get_eigenvals()[spin](kpoint, band) =
                    100.0 * spin + 10.0 * kpoint + band;
            }
        }
    }
    const auto previous = eigenvalue_snapshot(meanfield);
    meanfield.get_eigenvals()[1](1, 2) += 0.003;
    meanfield.get_eigenvals()[0](0, 0) -= 0.002;
    assert(std::abs(max_eigenvalue_change(meanfield, previous) - 0.003) <
           1.0e-14);
    assert(!qsgw_iteration_converged(4, 5, 1.0e-5, 1.0e-4));
    assert(qsgw_iteration_converged(5, 5, 1.0e-5, 1.0e-4));
}

void test_invalid_snapshot_or_tolerance_is_rejected()
{
    MeanField meanfield(1, 1, 1, 1, 1);
    auto snapshot = eigenvalue_snapshot(meanfield);
    snapshot.clear();
    assert_throws([&] { max_eigenvalue_change(meanfield, snapshot); });
    assert_throws([&] {
        qsgw_iteration_converged(1, 1, 0.0,
                                 std::numeric_limits<double>::quiet_NaN());
    });
}

} // namespace

int main()
{
    test_maximum_change_and_minimum_iteration_gate();
    test_invalid_snapshot_or_tolerance_is_rejected();
    std::cout << "test_qsgw_convergence: all tests passed\n";
    return 0;
}
