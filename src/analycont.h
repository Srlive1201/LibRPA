#pragma once
#include <vector>

#include "base_utility.h"

namespace LIBRPA
{

class AnalyContPade
{
private:
    int n_pars;
    std::vector<cplxdb> par_x;
    std::vector<cplxdb> par_y;

public:
    AnalyContPade(int n_pars_in,
                  const std::vector<cplxdb> &xs,
                  const std::vector<cplxdb> &data);
    cplxdb get(const cplxdb &x, const cplxdb &ref);
};

}
