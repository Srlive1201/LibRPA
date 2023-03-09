#include "../base_utility.h"
#include <cassert>

int main (int argc, char *argv[])
{

    assert(!is_complex<double>());
    assert(is_complex<std::complex<float>>());

    bool b;
    b = std::is_same<double, to_real<std::complex<double>>::type>::value;
    assert(b);

    return 0;
}
