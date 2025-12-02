#pragma once

#include "librpa.hpp"
#include "librpa_enums.h"

#include <string>

#include "../src/math/vector3_order.h"

namespace driver
{

// Runtime options specific to the driver
struct DriverParams
{
    std::string task;
    std::string input_dir = "./";

    // TODO: Move the following to the public LibrpaOptions class
    bool output_gw_spec_func;
    double sf_omega_start;
    double sf_omega_end;
    double sf_omega_step;
    double sf_gf_omega_shift;
    double sf_sigc_omega_shift;

    int sf_state_start;
    int sf_state_end;

    std::string format();

    DriverParams();
};

extern DriverParams driver_params;

extern const std::string input_filename;

// Types of each atom, read from structure file and also used to generate basis list
extern std::vector<int> atom_types;
extern size_t n_atoms;

// Dimension information, used across a few read_data functions
extern int n_spins;
extern int n_kpoints;
extern int n_ibz_kpoints;
extern int n_states;
extern int n_basis_wfc;

extern std::vector<std::pair<size_t, size_t>> local_atpair;
extern std::vector<librpa_int::Vector3_Order<double>> ibz_kpoints;

// Working handle
extern librpa::Handler h;

// Working runtime options
extern librpa::Options opts;

// TODO: consider move to public API
std::string format_runtime_options(const librpa::Options &opts) noexcept;

LibrpaTimeFreqGrid get_tfgrid_type(const std::string& grid_str);

std::string get_tfgrid_string(const LibrpaTimeFreqGrid& grid_type) noexcept;

inline LibrpaSwitch get_switch(bool switch_bool) noexcept
{
    return switch_bool ? LIBRPA_SWITCH_ON : LIBRPA_SWITCH_OFF;
}

inline bool get_bool(LibrpaSwitch switch_in) noexcept
{
    return switch_in == LIBRPA_SWITCH_ON;
}

LibrpaParallelRouting get_parallel_routing(const std::string& routing_str_low);

std::string get_routing_string(LibrpaParallelRouting routing);

}
