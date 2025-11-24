#include "driver_params.h"

#include "../src/io/global_io.h"

DriverParams::DriverParams():
    task("rpa"),
    input_dir(""),
    output_gw_spec_func(false),
    sf_omega_start(0.0),
    sf_omega_end(1.0),
    sf_omega_step(0.1),
    sf_gf_omega_shift(0.01),
    sf_sigc_omega_shift(0.01),
    sf_state_start(0),
    sf_state_end(10000)
{
}

void DriverParams::print()
{
    librpa_int::global::lib_printf("task = %s\n", task.c_str());
    librpa_int::global::lib_printf("input_dir = %s\n", input_dir.c_str());
    librpa_int::global::lib_printf("output_gw_spec_func = %L\n", output_gw_spec_func);
    if (output_gw_spec_func)
    {
        librpa_int::global::lib_printf("sf_omega_start = %f\n", sf_omega_start);
        librpa_int::global::lib_printf("sf_omega_end   = %f\n", sf_omega_end);
        librpa_int::global::lib_printf("sf_omega_step  = %f\n", sf_omega_step);
    }
}

DriverParams driver_params;
