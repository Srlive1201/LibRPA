#pragma once
#include <string>

struct DriverParams
{
    std::string input_dir;

    bool output_gw_spec_func;
    double sf_omega_start;
    double sf_omega_end;
    double sf_omega_step;
    double sf_gf_omega_shift;
    double sf_sigc_omega_shift;

    int sf_state_start;
    int sf_state_end;

    void print();

    DriverParams();
};

extern DriverParams driver_params;
