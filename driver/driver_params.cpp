#include "driver_params.h"

#include "utils_io.h"

DriverParams::DriverParams():
    input_dir(""),
    output_gw_spec_func(false),
    omega_sf_start(0.0),
    omega_sf_end(-1.0),
    omega_sf_step(0.1)
{
}

void DriverParams::print()
{
    LIBRPA::utils::lib_printf("input_dir = %s\n", input_dir.c_str());
    LIBRPA::utils::lib_printf("output_gw_spec_func = %L\n", output_gw_spec_func);
    if (output_gw_spec_func)
    {
        LIBRPA::utils::lib_printf("omega_sf_start = %f\n", omega_sf_start);
        LIBRPA::utils::lib_printf("omega_sf_end   = %f\n", omega_sf_end);
        LIBRPA::utils::lib_printf("omega_sf_step  = %f\n", omega_sf_step);
    }
}

DriverParams driver_params;
