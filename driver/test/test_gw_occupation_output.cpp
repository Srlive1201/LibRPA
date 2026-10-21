#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "librpa.h"

namespace
{
void require(const bool condition, const char *message)
{
    if (!condition) throw std::runtime_error(message);
}

void check_output(const int n_spins, const std::vector<double> &weights,
                  const std::vector<int> &mapping)
{
    const int nk = static_cast<int>(weights.size());
    const int n_states = 3;
    librpa::Handler handler(MPI_COMM_WORLD);
    handler.set_scf_dimension(n_spins, nk, n_states, n_states);
    const auto count = static_cast<std::size_t>(n_spins) * nk * n_states;
    std::vector<double> vxc(count, 0.0), vexx(count, 0.0), sigc(2 * count, 0.0);
    std::vector<double> energies(count), occupations_input(count);
    std::vector<double> expected_occupations, expected_energies;
    const double occupations[] = {n_spins == 1 ? 2.0 : 1.0, 0.4, 0.0};
    for (int spin = 0; spin < n_spins; ++spin)
    {
        for (int ik = 0; ik < nk; ++ik)
            for (int state = 0; state < n_states; ++state)
            {
                const auto index = (static_cast<std::size_t>(spin) * nk + ik) * n_states + state;
                // The public setter divides the supplied occupations by nk.
                occupations_input[index] = occupations[state] * weights[ik] * nk;
                energies[index] = 10.0 * spin + ik + 0.1 * state;
            }
    }
    handler.set_wg_ekb_efermi(n_spins, nk, n_states, occupations_input.data(), energies.data(),
                              0.0);
    const int n_output = mapping.empty() ? nk : static_cast<int>(mapping.size());
    std::vector<double> kfrac(3 * n_output, 0.0);
    for (int output = 0; output < n_output; ++output)
    {
        const int ik = mapping.empty() ? static_cast<int>(output) : mapping[output];
        for (int spin = 0; spin < n_spins; ++spin)
            for (int state = 0; state < n_states; ++state)
            {
                expected_occupations.push_back(occupations[state]);
                expected_energies.push_back(energies[(spin * nk + ik) * n_states + state]);
            }
    }
    LibrpaGwOutputInput view{};
    view.n_spins = n_spins;
    view.n_kpoints = nk;
    view.n_states = n_states;
    view.i_state_high = n_states;
    view.n_output_kpoints = n_output;
    view.fractional_k = kfrac.data();
    view.output_to_input_kpoint = mapping.empty() ? nullptr : mapping.data();
    view.vxc_ha = vxc.data();
    view.vexx_ha = vexx.data();
    view.sigc_ri = sigc.data();
    view.n_kpoint_weights = nk;
    view.kpoint_weights = weights.data();
    std::ofstream("energy_qp") << handler.format_gw_energy_table(view);
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
        {
            std::vector<double> loaded_energies(n_states), loaded_occupations(n_states);
            librpa_get_meanfield_k(handler.get_c_handler(), 0, spin, ik, loaded_energies.data(),
                                   loaded_occupations.data());
            for (int state = 0; state < n_states; ++state)
                require(loaded_occupations[state] ==
                            occupations_input[(spin * nk + ik) * n_states + state] * (1.0 / nk),
                        "output changed internal occupations");
        }
    for (const double weight : {0.0, -0.5, std::numeric_limits<double>::infinity(),
                                std::numeric_limits<double>::quiet_NaN()})
    {
        auto invalid_weights = weights;
        invalid_weights.front() = weight;
        view.kpoint_weights = invalid_weights.data();
        bool rejected = false;
        try
        {
            handler.format_gw_energy_table(view);
        }
        catch (const std::invalid_argument &)
        {
            rejected = true;
        }
        require(rejected, "invalid k-point weight was accepted");
    }
    view.kpoint_weights = weights.data();
    bool rejected = false;
    try
    {
        view.n_kpoint_weights = 0;
        handler.format_gw_energy_table(view);
    }
    catch (const std::invalid_argument &)
    {
        rejected = true;
    }
    require(rejected, "mismatched k-point weight count was accepted");
}
}  // namespace

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    librpa::init_global();
    const auto original_directory = std::filesystem::current_path();
    const auto test_directory = original_directory / "gw_occupation_output_test";
    bool directory_created = false;
    try
    {
        directory_created = std::filesystem::create_directory(test_directory);
        require(directory_created, "test output directory already exists");
        std::filesystem::current_path(test_directory);
        check_output(1, {0.125, 0.5, 0.375}, {0, 1, 1, 1, 1, 2, 2, 2});
        check_output(1, {0.125, 0.5, 0.375}, {});
        check_output(2, {0.25, 0.75}, {0, 1, 1, 1});
        check_output(1, {0.25, 0.25, 0.25, 0.25}, {});
        std::filesystem::current_path(original_directory);
        std::filesystem::remove(test_directory / "energy_qp");
        std::filesystem::remove(test_directory);
        std::cout << "GW occupation output checks passed\n";
        librpa::finalize_global();
        MPI_Finalize();
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
        librpa::finalize_global();
        MPI_Finalize();
        return 1;
    }
}
