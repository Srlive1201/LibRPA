#pragma once
#include <string>

struct DriverParams
{
    std::string input_dir;

    bool output_gw_spec_func;
    double omega_sf_start;
    double omega_sf_end;
    double omega_sf_step;

    void print();

    DriverParams();
};

extern DriverParams driver_params;
