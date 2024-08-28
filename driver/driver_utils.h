#pragma once
#include <vector>

std::vector<double> interpolate_dielec_func(int option, const std::vector<double> &frequencies_in,
                                            const std::vector<double> &df_in,
                                            const std::vector<double> &frequencies_target);
