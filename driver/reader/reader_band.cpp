#include <algorithm>
#include <cassert>
#include <fstream>
#include <iomanip>
#include <sstream>

#include "../../src/io/fs.h"
#include "../../src/io/global_io.h"
#include "../../src/utils/error.h"
#include "read_data.h"
namespace librpa::reader
{
using namespace librpa_int;
using namespace librpa_int::global;
using std::ifstream;
using std::stod;
using std::string;
void read_band_eigenvalues(ReaderContext &ctx, const string &dir_path)
{
    const int n_spins = ctx.state.n_spins;
    const int n_kpoints_band = ctx.state.n_kpoints_band;
    const int n_states = ctx.state.n_states;
    std::vector<double> eskb(n_spins * n_kpoints_band * n_states);
    std::vector<double> wskb(n_spins * n_kpoints_band * n_states);

    const int n_kb = n_kpoints_band * n_states;
    std::string s1, s2, s3, s4, s5;

    // Load occupation weights and eigenvalues
    for (int ik = 0; ik < n_kpoints_band; ik++)
    {
        std::stringstream ss;
        ss << dir_path << "band_KS_eigenvalue_k_" << std::setfill('0') << std::setw(5)
           << ik + 1 << ".txt";
        const auto file_path = ss.str();
        librpa_int::require_readable_file(file_path);
        ofs_myid << "Loading band eigenvalues from " << file_path << std::endl;
        ifstream infile;
        infile.open(file_path);
        for (int i_spin = 0; i_spin < n_spins; i_spin++)
        {
            for (int i_state = 0; i_state < n_states; i_state++)
            {
                infile >> s1 >> s2 >> s3 >> s4 >> s5;
                const int index = i_spin * n_kb + ik * n_states + i_state;
                wskb[index] = stod(s3);
                eskb[index] = stod(s4);
            }
        }
        infile.close();
    }
    ctx.h.set_band_occ_eigval(n_spins, n_kpoints_band, n_states, wskb.data(), eskb.data());
}

void read_band_meanfield_data(ReaderContext &ctx, const string &dir_path)
{
    using namespace librpa_int;
    using namespace librpa_int::global;
    auto &iks_band_eigvec_this = ctx.state.iks_band_eigvec_this;
    auto &n_basis_ao = ctx.state.n_basis_ao;
    auto &n_basis_wfc = ctx.state.n_basis_wfc;
    auto &n_kpoints_band = ctx.state.n_kpoints_band;
    auto &n_spins = ctx.state.n_spins;
    auto &n_states = ctx.state.n_states;
    using std::endl;

    if (ctx.state.n_kpoints_band == 0)
        throw LIBRPA_RUNTIME_ERROR(
            "Number of band k-points not set, run read_band_kpath_info first");

    iks_band_eigvec_this.clear();
    for (int ik = 0; ik < n_kpoints_band; ++ik)
        if (!reader_switch(ctx.opts.use_kpara_scf_eigvec) || ik % ctx.comm.nprocs == ctx.comm.myid)
            iks_band_eigvec_this.push_back(ik);

    read_band_eigenvalues(ctx, dir_path);

    // Load eigenvectors
    for (int ik = 0; ik < n_kpoints_band; ik++)
    {
        bool skip_this_ik = false;
        if (reader_switch(ctx.opts.use_kpara_scf_eigvec))
        {
            const auto it =
                std::find(iks_band_eigvec_this.cbegin(), iks_band_eigvec_this.cend(), ik);
            skip_this_ik = (it == iks_band_eigvec_this.cend());
        }
        if (skip_this_ik) continue;

        std::stringstream ss;
        ss << dir_path << "band_KS_eigenvector_k_" << std::setfill('0') << std::setw(5)
           << ik + 1 << ".txt";
        const auto file_path = ss.str();
        librpa_int::require_readable_file(file_path);

        ifstream infile;
        infile.open(file_path, std::ios::in | std::ios::binary);
        if (!infile.good())
            throw LIBRPA_RUNTIME_ERROR("Fail to open band eigenvector file " + file_path);
        else
            ofs_myid << "Loading band eigenvector file " + file_path << endl;

        size_t total_complex_comp =
            static_cast<size_t>(n_states) * static_cast<size_t>(n_basis_ao);  // for one component
        size_t total_complex_spin =
            static_cast<size_t>(n_states) * static_cast<size_t>(n_basis_wfc);
        size_t total_complex = total_complex_spin * n_spins;
        size_t bytes_doubles = total_complex * 2 * sizeof(double);

        std::vector<std::complex<double>> vecs(total_complex);
        infile.read(reinterpret_cast<char *>(vecs.data()), bytes_doubles);
        if (!infile || infile.gcount() != static_cast<ptrdiff_t>(bytes_doubles))
        {
            throw LIBRPA_RUNTIME_ERROR("Error: failed to read " + file_path);
        }

        const bool use_spinor_wfc = ctx.params.use_spinor_wfc;
        int n_soc = use_spinor_wfc ? 2 : 1;
        assert(n_soc * n_spins <= 2);
        for (int i_spin = 0; i_spin < n_spins; ++i_spin)
        {
            std::vector<std::complex<double>> vecs_sp(total_complex_spin);
            for (int ib = 0; ib < n_states; ++ib)
            {
                for (int iw = 0; iw < n_basis_ao; ++iw)
                {
                    for (int i_soc = 0; i_soc < n_soc; ++i_soc)
                    {
                        const size_t index_dst = i_soc * total_complex_comp + ib * n_basis_ao + iw;
                        size_t index_src;
                        if (use_spinor_wfc)
                        {
                            // NOTE: i_spin should be 0 for spinor-form wavefunction
                            assert(i_spin < 1);
                            index_src = ib * n_basis_ao * n_soc + iw * n_soc + i_soc;
                        }
                        else
                        {
                            index_src = i_spin * n_basis_ao * n_states + ib * n_basis_ao + iw;
                        }
                        vecs_sp[index_dst] = vecs[index_src];
                    }
                }
            }
            if (use_spinor_wfc)
            {
                ctx.h.set_wfc_band_spinor_packed(ik, n_states, n_basis_ao, vecs_sp.data(),
                                                 vecs_sp.data() + total_complex_comp);
            }
            else
            {
                ctx.h.set_wfc_band_packed(i_spin, ik, n_states, n_basis_ao, vecs_sp.data());
            }
        }

        infile.close();
    }
}

}  // namespace librpa::reader
