#include "../analycont.h"

#include <cassert>
#include "testutils.h"
#include <iostream>


void test_analycont_pade()
{
    using std::vector;
    using std::cout;

    int n = 6;
    // original function: y = -1 / (x - x0), x0 = {1.0, 1.0}
    const cplxdb x0 = {2.0, 2.0};
    std::vector<cplxdb> xs(n);
    std::vector<cplxdb> data(n);
    for (int i = 0; i < n; i++)
    {
        xs[i] = {0.0, static_cast<double>(i+1)};
        data[i] = -1.0 / (xs[i] - x0);
    }
    cout << "Input x: " << xs << "\n";
    cout << "Input y: " << data << "\n";

    LIBRPA::AnalyContPade pade(n, xs, data);
    // std::cout << data << "\n";
    const cplxdb test_x = {1.0, 1.0};
    const cplxdb ref = -1.0 / (test_x - x0);
    cplxdb test_y = pade.get(test_x);
    cout << "Reference: "<< ref << "\n";
    cout << "     Test: "<< test_y << "\n";
    assert(fequal(ref, test_y));
}

int main(int argc, char *argv[])
{
    test_analycont_pade();
}
