// Public API headers
#include "librpa_file_reader.hpp"

// Internal headers
#include "instance_manager.h"
#include "../utils/error.h"
#include "../../driver/driver.h"
#include "../../driver/read_data.h"
#include "../../driver/reader_coulomb.h"
#include "../../driver/reader_lri.h"

// Standard headers
#include <exception>
#include <filesystem>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace librpa
{
namespace
{

std::string normalize_input_dir(const std::string &path)
{
    if (path.empty())
    {
        throw LIBRPA_RUNTIME_ERROR("Input directory must not be empty");
    }
    return path.back() == '/' ? path : path + '/';
}

void require_input_file(const fs::path &path)
{
    if (!fs::is_regular_file(path))
    {
        throw LIBRPA_RUNTIME_ERROR("Required LibRPA input file is missing: " + path.string());
    }
}

template <typename Function>
void run_reader(const char *name, Function &&function)
{
    try
    {
        function();
    }
    catch (const std::exception &error)
    {
        throw LIBRPA_RUNTIME_ERROR(std::string("LibRPA file reader '") + name
                                   + "' failed: " + error.what());
    }
}

void read_scf_meanfield(const fs::path &path)
{
    librpa_int::MeanField meanfield;
    read_scf_occ_eigenvalues(path.string(), meanfield, false);

    driver::n_spins = meanfield.get_n_spins();
    driver::n_kpoints = meanfield.get_n_kpoints();
    driver::n_ibz_kpoints = driver::n_kpoints;
    driver::n_states = meanfield.get_n_states();
    driver::n_basis_wfc = meanfield.get_n_aos();
    driver::n_basis_ao = meanfield.get_n_aos();
    driver::n_spinor = meanfield.get_n_spinor();
    driver::iks_eigvec_this.resize(static_cast<std::size_t>(driver::n_kpoints));
    for (int ik = 0; ik != driver::n_kpoints; ++ik)
    {
        driver::iks_eigvec_this[static_cast<std::size_t>(ik)] = ik;
    }

    const std::size_t block_size = static_cast<std::size_t>(driver::n_kpoints)
                                   * static_cast<std::size_t>(driver::n_states);
    std::vector<double> eigenvalues(static_cast<std::size_t>(driver::n_spins) * block_size);
    std::vector<double> occupations(eigenvalues.size());
    for (int is = 0; is != driver::n_spins; ++is)
    {
        const auto &eigenvalues_spin
            = meanfield.get_eigenvals()[static_cast<std::size_t>(is)];
        const auto &occupations_spin
            = meanfield.get_weight()[static_cast<std::size_t>(is)];
        for (std::size_t index = 0; index != block_size; ++index)
        {
            eigenvalues[static_cast<std::size_t>(is) * block_size + index]
                = eigenvalues_spin.c[index];
            // The file reader normalizes weights, while Handler performs this
            // normalization when setting the dataset. Restore the file values.
            occupations[static_cast<std::size_t>(is) * block_size + index]
                = occupations_spin.c[index] * driver::n_kpoints;
        }
    }

    driver::h.set_scf_dimension(driver::n_spins, driver::n_kpoints,
                                driver::n_states, driver::n_basis_ao,
                                driver::n_spinor);
    driver::h.set_wg_ekb_efermi(driver::n_spins, driver::n_kpoints,
                               driver::n_states, occupations.data(),
                               eigenvalues.data(), meanfield.get_efermi());
}

} // namespace

std::shared_ptr<librpa_int::Dataset> read_dataset_from_files(
    MPI_Comm comm, const FileReaderOptions &opts)
{
    const std::string input_dir = normalize_input_dir(opts.input_dir);
    const fs::path input_path(input_dir);
    require_input_file(input_path / "band_out");
    require_input_file(input_path / "stru_out");

    driver::h.init(comm);
    driver::driver_params = driver::DriverParams{};
    driver::opts = Options{};
    driver::driver_params.input_dir = input_dir;
    driver::driver_params.cs_threshold = opts.cs_threshold;
    driver::driver_params.prefix_coul_cut = "coulomb_unshrinked_cut";
    driver::opts.parallel_routing = LIBRPA_ROUTING_LIBRI;
    driver::opts.vq_threshold = opts.coulomb_threshold;
    driver::opts.use_kpara_scf_eigvec = LIBRPA_SWITCH_OFF;

    read_scf_meanfield(input_path / "band_out");
    read_stru((input_path / "stru_out").string());
    if (fs::is_regular_file(input_path / "bz_sampling_out"))
    {
        read_bz_sampling((input_path / "bz_sampling_out").string());
    }
    else
    {
        read_bz_sampling_from_stru((input_path / "stru_out").string());
    }
    read_basis_wfc_aux(input_dir,
                       driver::driver_params.fn_basis,
                       driver::driver_params.fn_basis_wfc,
                       driver::driver_params.fn_basis_aux);

    auto pds = librpa_int::api::get_dataset_instance(driver::h.get_c_handler());
    const auto aux_sizes = pds->basis_aux.get_atom_nbs();
    if (aux_sizes.size() != pds->atoms.size())
    {
        throw LIBRPA_RUNTIME_ERROR(
            "LibRPA auxiliary-basis atom partition is inconsistent: "
            + std::to_string(aux_sizes.size()) + " basis blocks for "
            + std::to_string(pds->atoms.size()) + " atoms");
    }

    if (opts.read_ri)
    {
        const auto all_pairs
            = librpa_int::generate_atom_pair_from_nat(driver::n_atoms, false);
        run_reader("coulomb_unshrinked_cut", [&]
        {
            read_Vq_row(input_dir, "coulomb_unshrinked_cut",
                        opts.coulomb_threshold, all_pairs, true,
                        driver::driver_params.version_coul_reader, false);
        });
        run_reader("Cs_data", [&]
        {
            read_Cs(input_dir, opts.cs_threshold, all_pairs,
                    driver::driver_params.prefix_lri_coeff,
                    driver::driver_params.version_lri_reader);
        });
    }

    run_reader("coarse KS eigenvectors", [&]
    {
        if (read_eigenvector(input_dir, pds->mf, false) != 0)
        {
            throw LIBRPA_RUNTIME_ERROR("Failed to read coarse-grid KS eigenvectors");
        }
    });
    require_input_file(input_path / "velocity_matrix");
    run_reader("velocity_matrix", [&]
    {
        read_velocity((input_path / "velocity_matrix").string(), pds->mf,
                      pds->velocity_matrix);
    });
    if (opts.read_band_data)
    {
        require_input_file(input_path / "band_kpath_info");
        run_reader("band_kpath_info", [&]
        {
            read_band_kpath_info((input_path / "band_kpath_info").string());
        });
        run_reader("band meanfield", [&]
        {
            read_band_meanfield_data(input_dir);
        });
    }
    else
    {
        pds->kfrac_band_list = pds->pbc.kfrac_list;
        pds->mf_band = pds->mf;
    }

    return pds;
}

} // namespace librpa
