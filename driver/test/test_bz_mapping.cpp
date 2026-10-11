#include <cassert>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

#include "../../src/api/instance_manager.h"
#include "../../src/io/global_io.h"
#include "../../src/mpi/global_mpi.h"
#include "../../src/utils/constants.h"
#include "../reader/read_data.h"
#include "../reader/reader_context.h"

namespace
{
void write_sampling(const std::filesystem::path &path, const std::string &tail)
{
    std::ofstream out(path);
    out << std::setprecision(17);
    out << "1 1 3\n2 2\n";
    out << "1 " << 1.0 / 3 << " 0 0 0 0 0 0 1 1\n";
    out << "2 " << 2.0 / 3 << " 0 0 " << 1.0 / 3 << " 0 0 " << librpa_int::TWO_PI / 3 << " 2 2\n";
    out << tail;
}

std::string valid_tail()
{
    std::ostringstream out;
    out << std::setprecision(17);
    out << "full_kmap 3\n1 0 0 0 1\n2 0 0 " << librpa_int::TWO_PI / 3 << " 2\n3 0 0 "
        << 2 * librpa_int::TWO_PI / 3 << " 2\n";
    return out.str();
}

void check_mapping(const std::filesystem::path &file, const std::vector<int> &sizes,
                   const MPI_Comm comm = MPI_COMM_SELF, const bool explicit_mapping = true)
{
    librpa::Handler h(comm);
    const double lattice[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
    const double reciprocal[9] = {librpa_int::TWO_PI, 0, 0, 0, librpa_int::TWO_PI, 0, 0, 0,
                                  librpa_int::TWO_PI};
    h.set_latvec_and_G(lattice, reciprocal);
    const int n_scf = static_cast<int>(sizes.size());
    h.set_scf_dimension(1, n_scf, 1, 1, 1);
    h.set_atoms({0}, {0.0, 0.0, 0.0});
    const int identity[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
    h.set_symmetry_operations(1, 1, identity);
    librpa::reader::ReaderState state;
    state.n_kpoints = n_scf;
    librpa::reader::ReaderContext ctx{h, state, {}, {}, librpa_int::MpiCommHandler(comm, true)};
    librpa::reader::read_bz_sampling(ctx, file.string());
    const auto ds = librpa_int::api::get_dataset_instance(h.get_c_handler());
    const auto &pbc = ds->pbc;
    assert(pbc.klist.size() == sizes.size());
    assert(pbc.klist_full.size() == 3);
    for (int ik = 0; ik < n_scf; ++ik)
    {
        const auto expected_members = explicit_mapping ? sizes[ik] : 1;
        assert(pbc.map_irk_ks.at(pbc.klist_coul[ik]).size() ==
               static_cast<std::size_t>(expected_members));
        assert(std::abs(pbc.weight_k[ik] - sizes[ik] / 3.0) < 1e-12);
        assert(std::abs(pbc.weight_q[ik] - sizes[ik] / 3.0) < 1e-12);
    }
}
void check_explicit_api_validation()
{
    using namespace librpa_int;
    librpa::Handler h(MPI_COMM_SELF);
    const double lattice[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
    const double reciprocal[9] = {TWO_PI, 0, 0, 0, TWO_PI, 0, 0, 0, TWO_PI};
    h.set_latvec_and_G(lattice, reciprocal);
    h.set_scf_dimension(1, 1, 1, 1, 1);
    const std::vector<double> scf = {0, 0, 0};
    const std::vector<double> full = {0, 0, 0, TWO_PI / 3, 0, 0, 2 * TWO_PI / 3, 0, 0};
    h.set_irreducible_kgrids_kvec(3, 1, 1, scf, full, {0, 0, 0});
    auto ds = api::get_dataset_instance(h.get_c_handler());
    const auto saved_full = ds->pbc.klist_full;
    const auto saved_weights = ds->pbc.weight_k;
    ds->symmetry_context.available = true;
    for (const auto &bad : std::vector<std::vector<double>>{
             {0, 0, 0, 0.1 * TWO_PI, 0, 0, 0.2 * TWO_PI, 0, 0},
             {0, 0, 0, TWO_PI, 0, 0, 2 * TWO_PI / 3, 0, 0},
             {TWO_PI / 6, 0, 0, TWO_PI / 2, 0, 0, 5 * TWO_PI / 6, 0, 0}})
    {
        bool rejected = false;
        try
        {
            h.set_irreducible_kgrids_kvec(3, 1, 1, scf, bad, {0, 0, 0});
        }
        catch (const std::exception &)
        {
            rejected = true;
        }
        assert(rejected);
        assert(ds->pbc.klist_full == saved_full);
        assert(ds->pbc.weight_k == saved_weights);
        assert(ds->symmetry_context.available);
    }
    h.set_irreducible_kgrids_kvec(3, 1, 1, scf, full, {0, 0, 0});
    assert(!ds->symmetry_context.available);
}

void check_shifted_nonorthogonal_mapping(const std::filesystem::path &file)
{
    using namespace librpa_int;
    librpa::Handler h(MPI_COMM_SELF);
    const double lattice[9] = {2, 0, 0, 1, 3, 0, 0, 0, 4};
    const double reciprocal[9] = {TWO_PI / 2, -TWO_PI / 6, 0, 0, TWO_PI / 3, 0, 0, 0, TWO_PI / 4};
    h.set_latvec_and_G(lattice, reciprocal);
    h.set_scf_dimension(1, 2, 1, 1, 1);
    const auto k = [](const double x)
    {
        return Vector3_Order<double>{x * TWO_PI / 2, -x * TWO_PI / 6 + 0.2 * TWO_PI / 3,
                                     0.1 * TWO_PI / 4};
    };
    std::ofstream out(file);
    out << std::setprecision(17) << "3 1 1\n2 2\n";
    for (int ik = 0; ik < 2; ++ik)
    {
        const double x = 1.0 / 6 + ik / 3.0;
        const auto v = k(x);
        out << ik + 1 << ' ' << (ik + 1) / 3.0 << ' ' << x << " 0.2 0.1 " << v.x << ' ' << v.y
            << ' ' << v.z << ' ' << ik + 1 << ' ' << ik + 1 << '\n';
    }
    out << "full_kmap 3\n";
    for (int ik = 0; ik < 3; ++ik)
    {
        const auto v = k(1.0 / 6 + ik / 3.0);
        out << ik + 1 << ' ' << v.x << ' ' << v.y << ' ' << v.z << ' ' << (ik == 0 ? 1 : 2) << '\n';
    }
    out.close();
    librpa::reader::ReaderState state;
    state.n_kpoints = 2;
    librpa::reader::ReaderContext ctx{h, state, {}, {}, MpiCommHandler(MPI_COMM_SELF, true)};
    librpa::reader::read_bz_sampling(ctx, file.string());
    const auto ds = api::get_dataset_instance(h.get_c_handler());
    assert(ds->pbc.klist_full.size() == 3);
    assert(std::abs(ds->pbc.kfrac_list[0].x - 1.0 / 6) < 1e-12);
    assert(std::abs(ds->pbc.kfrac_list[1].x - 0.5) < 1e-12);
    assert(std::abs(ds->pbc.weight_k[0] - 1.0 / 3) < 1e-12);
}
}  // namespace

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    librpa_int::global::init_global_mpi(MPI_COMM_WORLD);
    librpa_int::global::init_global_io();
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    const auto stamp = std::chrono::steady_clock::now().time_since_epoch().count();
    const auto directory =
        std::filesystem::temp_directory_path() /
        ("librpa_bz_mapping_" + std::to_string(stamp) + "_" + std::to_string(rank));
    assert(std::filesystem::create_directory(directory));
    const auto file = directory / "bz_sample.txt";
    check_explicit_api_validation();
    // Reduced legacy files still use the original inference path.
    for (const std::string tail : std::vector<std::string>{
             "", "2\n1 0.3333333333333333\n2 0.6666666666666666\n", "legacy summary\n"})
    {
        write_sampling(file, tail);
        check_mapping(file, {1, 2}, MPI_COMM_SELF, false);
    }
    write_sampling(file, valid_tail());
    check_mapping(file, {1, 2});
    for (const std::string tail : std::vector<std::string>{
             "full_kmap 2\n", "full_kmap 3\n1 0 0 0 3\n",
             "full_kmap 3\n1 0 0 0 1\n2 0 0 0 2\n3 0 0 0 2\n", valid_tail() + "extra\n"})
    {
        write_sampling(file, tail);
        bool rejected = false;
        try
        {
            check_mapping(file, {1, 2});
        }
        catch (const std::exception &)
        {
            rejected = true;
        }
        assert(rejected);
    }
    // Reject a complete, nonduplicated mapping that assigns the wrong star.
    std::string wrong_star = valid_tail();
    wrong_star.replace(wrong_star.find("1 0 0 0 1"), 9, "1 0 0 0 2");
    write_sampling(file, wrong_star);
    bool rejected = false;
    try
    {
        check_mapping(file, {1, 2});
    }
    catch (const std::exception &)
    {
        rejected = true;
    }
    assert(rejected);
    // Folded weights must agree with the explicitly assigned star sizes.
    std::string wrong_weights = valid_tail();
    write_sampling(file, wrong_weights);
    {
        std::ifstream input(file);
        std::string content((std::istreambuf_iterator<char>(input)), {});
        const auto first_weight = content.find("0.33333333333333331");
        const auto second_weight = content.find("0.66666666666666663");
        assert(first_weight != std::string::npos && second_weight != std::string::npos);
        content.replace(second_weight, 19, "0.5");
        content.replace(first_weight, 19, "0.5");
        std::ofstream(file) << content;
    }
    rejected = false;
    try
    {
        check_mapping(file, {1, 2});
    }
    catch (const std::exception &)
    {
        rejected = true;
    }
    assert(rejected);
    // Reject an overflowing mesh before allocation.
    std::ofstream(file) << "2147483647 2 2\n";
    rejected = false;
    try
    {
        check_mapping(file, {1, 2});
    }
    catch (const std::exception &)
    {
        rejected = true;
    }
    assert(rejected);
    // Existing full-grid files require no extension and retain unit stars.
    {
        std::ofstream out(file);
        out << std::setprecision(17) << "1 1 3\n3 3\n";
        for (int ik = 0; ik < 3; ++ik)
            out << ik + 1 << ' ' << 1.0 / 3 << " 0 0 " << ik / 3.0 << " 0 0 "
                << ik * librpa_int::TWO_PI / 3 << ' ' << ik + 1 << ' ' << ik + 1 << '\n';
    }
    check_mapping(file, {1, 1, 1}, MPI_COMM_SELF, false);
    {
        std::ofstream out(file, std::ios::app);
        out << "3\n1 0.3333333333333333\n2 0.3333333333333333\n3 0.3333333333333333\n";
    }
    check_mapping(file, {1, 1, 1}, MPI_COMM_SELF, false);
    check_shifted_nonorthogonal_mapping(file);
    // A non-root file-open failure must reach every rank before any reader API collective.
    write_sampling(file, valid_tail());
    int size = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size > 1)
    {
        bool all_rejected = false;
        try
        {
            check_mapping(rank == 1 ? directory / "missing.txt" : file, {1, 2}, MPI_COMM_WORLD);
        }
        catch (const std::exception &)
        {
            all_rejected = true;
        }
        assert(all_rejected);
    }
    std::filesystem::remove(file);
    std::filesystem::remove(directory);
    librpa_int::global::finalize_global_io();
    librpa_int::global::finalize_global_mpi();
    MPI_Finalize();
}
