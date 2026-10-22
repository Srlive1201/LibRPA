#include <cassert>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "../../driver/tasks/task_helper.h"

int main()
{
    using namespace librpa_int;

    MeanField mf(1, 2, 1, 1);
    mf.get_weight()[0](0, 0) = 0.5;
    mf.get_weight()[0](1, 0) = 1.5;

    std::vector<Vector3_Order<double>> kfrac{
        {0.0, 0.0, 0.0}, {0.25, 0.0, 0.0}, {0.5, 0.0, 0.0}, {0.75, 0.0, 0.0}};
    std::vector<matrix> vxc(1);
    vxc[0].create(2, 1);
    vxc[0].zero_out();
    std::vector<double> vexx(2, 0.0);
    std::vector<cplxdb> sigc(2, 0.0);

    write_energy_qp(mf, kfrac, {0, 1, 1, 0}, vxc, vexx, sigc, 2, 0, 1,
                    {0.25, 0.75});

    std::ifstream ifs("energy_qp");
    assert(ifs.good());
    std::vector<double> occupations;
    std::string line;
    while (std::getline(ifs, line))
    {
        std::istringstream iss(line);
        int state = 0;
        double occupation = 0.0;
        double e_gs = 0.0;
        double e_qp = 0.0;
        if (iss >> state >> occupation >> e_gs >> e_qp) occupations.push_back(occupation);
    }
    std::remove("energy_qp");

    assert(occupations.size() == 4);
    for (const double occupation : occupations)
        assert(std::abs(occupation - 2.0) < 1.0e-12);
    return 0;
}
