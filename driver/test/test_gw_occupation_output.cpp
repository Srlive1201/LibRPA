#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "../tasks/task_helper.h"

namespace
{
void require(const bool condition, const char *message)
{
    if (!condition) throw std::runtime_error(message);
}

void check_output(const int n_spins, const std::vector<double> &weights,
                  const std::vector<int> &mapping)
{
    using namespace librpa_int;
    const int nk = static_cast<int>(weights.size());
    const int n_states = 3;
    MeanField mf(n_spins, nk, n_states, n_states);
    std::vector<matrix> vxc(n_spins);
    std::vector<double> vexx(n_spins * nk * n_states, 0.0);
    std::vector<cplxdb> sigc(n_spins * nk * n_states, 0.0);
    std::vector<double> expected_occupations, expected_energies;
    const double occupations[] = {n_spins == 1 ? 2.0 : 1.0, 0.4, 0.0};
    for (int spin = 0; spin < n_spins; ++spin)
    {
        vxc[spin].create(nk, n_states);
        vxc[spin].zero_out();
        for (int ik = 0; ik < nk; ++ik)
            for (int state = 0; state < n_states; ++state)
            {
                mf.get_weight()[spin](ik, state) = occupations[state] * weights[ik];
                mf.get_eigenvals()[spin](ik, state) = 10.0 * spin + ik + 0.1 * state;
            }
    }
    std::vector<Vector3_Order<double>> kfrac(mapping.empty() ? nk : mapping.size());
    for (std::size_t output = 0; output < kfrac.size(); ++output)
    {
        const int ik = mapping.empty() ? static_cast<int>(output) : mapping[output];
        for (int spin = 0; spin < n_spins; ++spin)
            for (int state = 0; state < n_states; ++state)
            {
                expected_occupations.push_back(occupations[state]);
                expected_energies.push_back(mf.get_eigenvals()[spin](ik, state));
            }
    }
    write_energy_qp(mf, kfrac, mapping, vxc, vexx, sigc, nk, 0, n_states, weights);
    std::ifstream input("energy_qp");
    require(input.good(), "energy_qp was not written");
    std::size_t row = 0;
    std::string line;
    while (std::getline(input, line))
    {
        std::istringstream fields(line);
        int state;
        double occupation, energy, qp_energy;
        if (!(fields >> state >> occupation >> energy >> qp_energy)) continue;
        require(row < expected_occupations.size(), "unexpected output row");
        require(state == static_cast<int>(row % n_states) + 1, "state index changed");
        require(std::abs(occupation - expected_occupations[row]) < 1.0e-12,
                "occupation was not restored using its mapped k-point weight");
        require(std::abs(energy - expected_energies[row]) < 1.0e-10, "mean-field energy changed");
        require(std::abs(qp_energy - expected_energies[row]) < 1.0e-10, "QP energy changed");
        ++row;
    }
    require(row == expected_occupations.size(), "missing output rows");
    for (int spin = 0; spin < n_spins; ++spin)
        for (int ik = 0; ik < nk; ++ik)
            for (int state = 0; state < n_states; ++state)
                require(mf.get_weight()[spin](ik, state) == occupations[state] * weights[ik],
                        "output changed internal occupations");
    bool rejected = false;
    try
    {
        write_energy_qp(mf, kfrac, mapping, vxc, vexx, sigc, nk, 0, n_states, {});
    }
    catch (const std::invalid_argument &)
    {
        rejected = true;
    }
    require(rejected, "mismatched k-point weight count was accepted");
}
}  // namespace

int main()
{
    const auto original_directory = std::filesystem::current_path();
    const auto test_directory = original_directory / "gw_occupation_output_test";
    bool directory_created = false;
    try
    {
        directory_created = std::filesystem::create_directory(test_directory);
        require(directory_created, "test output directory already exists");
        std::filesystem::current_path(test_directory);
        require(occupation_number_from_kpoint_weight(0.25, 0.125) == 2.0,
                "single-star occupation conversion failed");
        for (const double weight : {0.0, -0.5, std::numeric_limits<double>::infinity(),
                                    std::numeric_limits<double>::quiet_NaN()})
        {
            bool rejected = false;
            try
            {
                occupation_number_from_kpoint_weight(1.0, weight);
            }
            catch (const std::invalid_argument &)
            {
                rejected = true;
            }
            require(rejected, "invalid k-point weight was accepted");
        }
        check_output(1, {0.125, 0.5, 0.375}, {0, 1, 1, 1, 1, 2, 2, 2});
        check_output(1, {0.125, 0.5, 0.375}, {});
        check_output(2, {0.25, 0.75}, {0, 1, 1, 1});
        check_output(1, {0.25, 0.25, 0.25, 0.25}, {});
        std::filesystem::current_path(original_directory);
        std::filesystem::remove(test_directory / "energy_qp");
        std::filesystem::remove(test_directory);
        std::cout << "GW occupation output checks passed\n";
        return 0;
    }
    catch (const std::exception &error)
    {
        std::filesystem::current_path(original_directory);
        if (directory_created)
        {
            std::error_code ignored;
            std::filesystem::remove(test_directory / "energy_qp", ignored);
            std::filesystem::remove(test_directory, ignored);
        }
        std::cerr << "GW occupation output test failed: " << error.what() << '\n';
        return 1;
    }
}
