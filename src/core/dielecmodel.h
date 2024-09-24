#pragma once
#include <functional>
#include <vector>

namespace librpa_int {

//! double-dispersion Havriliak-Negami model
struct DoubleHavriliakNegami
{
    static const int d_npar;
    static const std::function<double(double, const std::vector<double> &)> func_imfreq;
    static const std::function<void(std::vector<double> &, double, const std::vector<double> &)>
        grad_imfreq;
};

std::vector<double> interpolate_dielec_func(int option, const std::vector<double> &frequencies_in,
                                            const std::vector<double> &df_in,
                                            const std::vector<double> &frequencies_target);
}

#include "../meanfield.h"

class diele_func
{
   private:
    // ( alpha, beta, omega )
    std::vector<std::vector<std::vector<std::complex<double>>>> head;
    std::vector<std::vector<std::vector<std::complex<double>>>> wing;
    const MeanField &meanfield_df;
    const std::vector<double> &omega;

   public:
    diele_func(const MeanField &mf, const std::vector<double> &frequencies_target)
        : meanfield_df(mf), omega(frequencies_target) {};
    ~diele_func() {};
    // All calculation in unit: Ang and eV.
    void cal_head();
    void cal_wing();
    double cal_factor();
    void init_head();
    void test_head();
};
