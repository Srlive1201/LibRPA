#include "../utils/base_utility.h"
#include <cassert>

using namespace librpa_int;

int main (int argc, char *argv[])
{

    assert(!is_complex<double>());
    assert(is_complex<std::complex<float>>());

    bool b;
    b = std::is_same<double, to_real<std::complex<double>>::type>::value;
    assert(b);

    // lexico comparison for pair
    auto row_fast = FastLess<std::pair<int, int>>{true}; // row index runs faster
    assert(row_fast({1, 0}, {0, 1}));
    auto col_fast = FastLess<std::pair<int, int>>{false};
    assert(col_fast({0, 1}, {1, 0}));

    return 0;
}
