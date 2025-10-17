#include "../analycont.h"
#include "../constants.h"
#include "../stl_io_helper.h"

#include <cassert>
#include "testutils.h"
#include <iostream>

using std::vector;
using std::cout;
using std::endl;

static void initialize_1_over_x_minus_x0_data(int n, const cplxdb &x0, vector<cplxdb> &xs,
                                              vector<cplxdb> &data)
{
    xs.resize(n);
    data.resize(n);
    for (int i = 0; i < n; i++)
    {
        xs[i] = {0.0, static_cast<double>(i+1)};
        data[i] = -1.0 / (xs[i] - x0);
    }
}

void test_analycont_pade()
{
    int n = 6;
    // original function: y = -1 / (x - x0), x0 = {2.0, 2.0}
    const cplxdb x0 = {2.0, 2.0};
    std::vector<cplxdb> xs(n);
    std::vector<cplxdb> data(n);
    initialize_1_over_x_minus_x0_data(n, x0, xs, data);
    cout << "Input x: " << xs << endl;
    cout << "Input y: " << data << endl;

    LIBRPA::AnalyContPade pade(n, xs, data);
    // std::cout << data << endl;
    const cplxdb test_x = {1.0, 1.0};
    const cplxdb ref = -1.0 / (test_x - x0);
    cplxdb test_y = pade.get(test_x);
    cout << "Reference: "<< ref << endl;
    cout << "     Test: "<< test_y << endl;
    assert(fequal(ref, test_y));
}

void test_get_spectfunc_acpade()
{
    int n = 6;
    // original function: y = -1 / (x - x0), x0 = {2.0, 2.0}
    const cplxdb x0 = {2.0, 2.0};
    std::vector<cplxdb> xs(n);
    std::vector<cplxdb> data(n);
    initialize_1_over_x_minus_x0_data(n, x0, xs, data);
    LIBRPA::AnalyContPade pade(n, xs, data);

    vector<cplxdb> omegas(2);
    omegas[0] = {0.0, 0.0};
    omegas[1] = {2.0, 3.0};

    // assume e_ks = v_xc = v_exx = 0, and then use the function above, we have
    // A(w) = - 1 / PI * Im [ w + 1 / (w - 2 - 2i) ]^{-1}
    // and
    // A(w=0) = - 1 / PI * Im [-2 - 2i]
    // A(w=2+3i) = - 1 / PI * Im [2 + 2i]^{-1}
    vector<double> ref_sf(2);
    ref_sf[0] = 2.0 / PI;
    ref_sf[1] = 0.25 / PI;

    auto sf = LIBRPA::get_specfunc(pade, omegas, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    assert(fequal_array(2, sf.data(), ref_sf.data(), false));

    // same results if we have v_xc cancels e_ks and v_exx
    sf = LIBRPA::get_specfunc(pade, omegas, 0.0, 1.0, 2.0, 1.0, 0.0, 0.0);
    assert(fequal_array(2, sf.data(), ref_sf.data(), false));

    // same results if we shift the omegas and e_ks by reference (2.0)
    omegas[0] = {2.0, -1.0};
    omegas[1] = {4.0, 2.0};
    sf = LIBRPA::get_specfunc(pade, omegas, 2.0, 1.0 + 2.0, 2.0, 1.0, 1.0, 1.0);
    assert(fequal_array(2, sf.data(), ref_sf.data(), false));
}

int main(int argc, char *argv[])
{
    test_analycont_pade();
    test_get_spectfunc_acpade();
}
