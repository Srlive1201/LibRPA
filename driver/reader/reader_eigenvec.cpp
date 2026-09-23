#include "reader_eigenvec.h"

#include <algorithm>
#include <cassert>
#include <complex>
#include <cstdint>
#include <fstream>
#include <ios>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "../../src/io/fs.h"
#include "../../src/io/global_io.h"
#include "../../src/core/meanfield_mpi.h"
#include "../../src/utils/profiler.h"
#include "reader_context.h"

namespace librpa::reader
{

namespace
{

constexpr std::int32_t READER_EIGENVEC_V1_MARKER = -12345679;
constexpr std::int32_t EIGENVEC_V1_KIND_COMPLEX_DOUBLE = 28;

static_assert(sizeof(std::complex<double>) == 2 * sizeof(double),
              "KS eigenvector v1 expects std::complex<double> as two doubles");

struct KpointBlock
{
    std::int32_t ik;
    std::int64_t offset;
};

struct WfcShape
{
    int nspin;
    int nsoc;
    int nband;
    int nao;
    int nkpoints;
    bool use_spinor;
};

int check_KS_file_version(const std::string &file_path)
{
    std::ifstream infile(file_path, std::ios::in | std::ios::binary);
    if (!infile.good()) return false;
    std::int32_t marker = 0;
    infile.read(reinterpret_cast<char *>(&marker), sizeof(marker));
    // Legacy text file
    if (!infile.good() || marker > 0) return 0;
    // Version 1 binary
    if (marker == READER_EIGENVEC_V1_MARKER) return 1;
    throw std::runtime_error("Unsupported KS eigenvector file marker: " + std::to_string(marker) +
                             " in " + file_path);
}

template <typename T>
bool read_one(std::ifstream &infile, T &value)
{
    infile.read(reinterpret_cast<char *>(&value), sizeof(value));
    return infile.good();
}

std::size_t wfc_index(const WfcShape &shape, const int is, const int isoc, const int ib,
                      const int iw)
{
    const auto nbao = static_cast<std::size_t>(shape.nband) * shape.nao;
    const auto n = static_cast<std::size_t>(shape.nsoc) * nbao;
    return static_cast<std::size_t>(is) * n + static_cast<std::size_t>(isoc) * nbao +
           static_cast<std::size_t>(ib) * shape.nao + static_cast<std::size_t>(iw);
}

bool selected_ik(const std::vector<int> *iks_selected, const int ik)
{
    return iks_selected == nullptr ||
           std::find(iks_selected->cbegin(), iks_selected->cend(), ik) != iks_selected->cend();
}

bool selected_input_ik(ReaderContext &ctx, const int ik)
{
    if (!reader_switch(ctx.opts.use_kpara_scf_eigvec)) return true;
    return std::find(ctx.state.iks_eigvec_this.cbegin(), ctx.state.iks_eigvec_this.cend(), ik) !=
           ctx.state.iks_eigvec_this.cend();
}

std::vector<std::string> eigenvector_files(ReaderContext &ctx, const std::string &dir_path)
{
    return librpa_int::discover_files_with_prefix(
        dir_path, ctx.params.prefix_eigvecs_scf);
}

template <typename ShouldReadIk, typename StoreIk>
int read_legacy_text_file(const std::string &file_path, const WfcShape &shape,
                          ShouldReadIk should_read_ik, StoreIk store_ik,
                          const LegacyTextWfcOrder text_order =
                              LegacyTextWfcOrder::BasisSpinorBandSpin)
{
    std::ifstream infile(file_path);
    if (!infile.good()) return 1;

    const auto n = static_cast<std::size_t>(shape.nspin) * shape.nsoc * shape.nband * shape.nao;
    std::vector<double> re(n);
    std::vector<double> im(n);

    std::string rvalue, ivalue, kstr;
    while (infile >> kstr)
    {
        const int ik = std::stoi(kstr) - 1;
        const bool keep_ik = should_read_ik(ik);
        std::fill(re.begin(), re.end(), 0.0);
        std::fill(im.begin(), im.end(), 0.0);

        if (text_order == LegacyTextWfcOrder::BasisSpinorBandSpin)
        {
            for (int iw = 0; iw != shape.nao; ++iw)
            {
                for (int isoc = 0; isoc != shape.nsoc; ++isoc)
                {
                    for (int ib = 0; ib != shape.nband; ++ib)
                    {
                        for (int is = 0; is != shape.nspin; ++is)
                        {
                            if (!(infile >> rvalue >> ivalue)) return 1;
                            if (!keep_ik) continue;
                            const auto dst = wfc_index(shape, is, isoc, ib, iw);
                            re[dst] = std::stod(rvalue);
                            im[dst] = std::stod(ivalue);
                        }
                    }
                }
            }
        }
        else
        {
            for (int is = 0; is != shape.nspin; ++is)
            {
                for (int iwfc = 0; iwfc != shape.nao * shape.nsoc; ++iwfc)
                {
                    const int isoc = iwfc % shape.nsoc;
                    const int iw = iwfc / shape.nsoc;
                    for (int ib = 0; ib != shape.nband; ++ib)
                    {
                        if (!(infile >> rvalue >> ivalue)) return 1;
                        if (!keep_ik) continue;
                        const auto dst = wfc_index(shape, is, isoc, ib, iw);
                        re[dst] = std::stod(rvalue);
                        im[dst] = std::stod(ivalue);
                    }
                }
            }
        }
        if (keep_ik) store_ik(ik, re, im);
    }
    return 0;
}

template <typename ShouldReadIk, typename StoreIk>
int read_binary_v1_file(const std::string &file_path, const WfcShape &shape,
                        ShouldReadIk should_read_ik, StoreIk store_ik)
{
    std::ifstream infile(file_path, std::ios::in | std::ios::binary);
    if (!infile.good()) return 1;

    std::int32_t marker = 0;
    std::int32_t kind_raw = 0;
    std::int32_t nkpoints_file = 0;
    std::int32_t nspins_file = 0;
    std::int32_t nstates_file = 0;
    std::int32_t nbasis_wfc_file = 0;
    if (!read_one(infile, marker) || !read_one(infile, kind_raw) ||
        !read_one(infile, nkpoints_file) || !read_one(infile, nspins_file) ||
        !read_one(infile, nstates_file) || !read_one(infile, nbasis_wfc_file))
    {
        return 1;
    }

    const int nbasis_wfc = shape.nao * shape.nsoc;
    if (marker != READER_EIGENVEC_V1_MARKER || kind_raw != EIGENVEC_V1_KIND_COMPLEX_DOUBLE ||
        nkpoints_file < 0 || nspins_file != shape.nspin || nstates_file != shape.nband ||
        nbasis_wfc_file != nbasis_wfc)
    {
        return 1;
    }

    std::vector<KpointBlock> blocks(static_cast<std::size_t>(nkpoints_file));
    for (auto &block : blocks)
    {
        if (!read_one(infile, block.ik) || !read_one(infile, block.offset))
        {
            return 1;
        }
        if (block.offset < 0) return 1;
    }

    const auto n = static_cast<std::size_t>(shape.nspin) * shape.nsoc * shape.nband * shape.nao;
    std::vector<std::complex<double>> block_data(n);

    for (const auto &block : blocks)
    {
        const int ik = block.ik - 1;
        if (ik < 0 || ik >= shape.nkpoints) return 1;
        if (!should_read_ik(ik)) continue;

        infile.clear();
        infile.seekg(static_cast<std::streamoff>(block.offset), std::ios::beg);
        if (!infile.good()) return 1;

        infile.read(reinterpret_cast<char *>(block_data.data()),
                    static_cast<std::streamsize>(block_data.size() * sizeof(std::complex<double>)));
        if (!infile.good()) return 1;
        store_ik(ik, block_data);
    }
    return 0;
}

void set_wfc_packed(ReaderContext &ctx, const WfcShape &shape, const int ik,
                           const std::vector<std::complex<double>> &block_data)
{
    const int nbao = shape.nband * shape.nao;
    for (int is = 0; is != shape.nspin; ++is)
    {
        const auto *spin_block =
            block_data.data() + static_cast<std::size_t>(is) * shape.nsoc * nbao;
        if (shape.use_spinor)
        {
            assert(is == 0);
            ctx.h.set_wfc_spinor_packed(ik, shape.nband, shape.nao, spin_block,
                                            spin_block + nbao);
        }
        else
        {
            ctx.h.set_wfc_packed(is, ik, shape.nband, shape.nao, spin_block);
        }
    }
}

void set_wfc(ReaderContext &ctx, const WfcShape &shape, const int ik, const std::vector<double> &re,
                    const std::vector<double> &im)
{
    const int nbao = shape.nband * shape.nao;
    const int n = shape.nsoc * nbao;
    for (int is = 0; is != shape.nspin; ++is)
    {
        if (shape.use_spinor)
        {
            assert(is == 0);
            ctx.h.set_wfc_spinor(ik, shape.nband, shape.nao, re.data(), im.data(),
                                     re.data() + nbao, im.data() + nbao);
        }
        else
        {
            ctx.h.set_wfc(is, ik, shape.nband, shape.nao,
                              re.data() + static_cast<std::size_t>(is) * n,
                              im.data() + static_cast<std::size_t>(is) * n);
        }
    }
}

void set_meanfield_wfc_packed(librpa_int::MeanField &mf, const WfcShape &shape, const int ik,
                              const std::vector<std::complex<double>> &block_data)
{
    const int nbao = shape.nband * shape.nao;
    for (int is = 0; is != shape.nspin; ++is)
    {
        const auto *spin_block =
            block_data.data() + static_cast<std::size_t>(is) * shape.nsoc * nbao;
        if (shape.use_spinor)
        {
            assert(is == 0);
            auto &wfcs_up = mf.get_eigenvectors()[0][0][ik];
            auto &wfcs_dn = mf.get_eigenvectors()[0][1][ik];
            wfcs_up.create(shape.nband, shape.nao);
            wfcs_dn.create(shape.nband, shape.nao);
            std::copy(spin_block, spin_block + nbao, wfcs_up.c);
            std::copy(spin_block + nbao, spin_block + 2 * nbao, wfcs_dn.c);
        }
        else
        {
            auto &wfc = mf.get_eigenvectors()[is][0][ik];
            wfc.create(shape.nband, shape.nao);
            std::copy(spin_block, spin_block + nbao, wfc.c);
        }
    }
}

void set_meanfield_wfc(librpa_int::MeanField &mf, const WfcShape &shape, const int ik,
                       const std::vector<double> &re, const std::vector<double> &im)
{
    const int nbao = shape.nband * shape.nao;
    const int n = shape.nsoc * nbao;
    for (int is = 0; is != shape.nspin; ++is)
    {
        if (shape.use_spinor)
        {
            assert(is == 0);
            auto &wfcs_up = mf.get_eigenvectors()[0][0][ik];
            auto &wfcs_dn = mf.get_eigenvectors()[0][1][ik];
            wfcs_up.create(shape.nband, shape.nao);
            wfcs_dn.create(shape.nband, shape.nao);
            for (int i = 0; i != nbao; ++i)
            {
                wfcs_up.c[i] = std::complex<double>(re[i], im[i]);
                wfcs_dn.c[i] = std::complex<double>(re[nbao + i], im[nbao + i]);
            }
        }
        else
        {
            auto &wfc = mf.get_eigenvectors()[is][0][ik];
            wfc.create(shape.nband, shape.nao);
            for (int i = 0; i != nbao; ++i)
            {
                const auto src = static_cast<std::size_t>(is) * n + i;
                wfc.c[i] = std::complex<double>(re[src], im[src]);
            }
        }
    }
}

void validate_kblacs_wfc_layout(const librpa_int::MeanField &mf,
                                const librpa_int::KPointBlacsParallelContext &kblacs_ctxt,
                                const librpa_int::ArrayDesc &desc_wfc)
{
    if (!kblacs_ctxt.is_initialized() || !desc_wfc.is_initialized())
        throw std::runtime_error("k-BLACS eigenvector reader requires initialized context and descriptor");
    if (kblacs_ctxt.n_kpoints() != mf.get_n_kpoints() ||
        desc_wfc.m() != mf.get_n_aos() || desc_wfc.n() != mf.get_n_states() ||
        desc_wfc.ictxt() != kblacs_ctxt.blacs_h.ictxt)
        throw std::runtime_error("k-BLACS eigenvector reader got inconsistent dimensions");
    const int block_ao = librpa_int::get_capped_blacs_block_size(
        mf.get_n_aos(), librpa_int::wfc_gemm_block_size_opt, kblacs_ctxt.blacs_h);
    const int block_band = librpa_int::get_capped_blacs_block_size(
        mf.get_n_states(), librpa_int::wfc_gemm_block_size_opt, kblacs_ctxt.blacs_h);
    if (desc_wfc.mb() != block_ao || desc_wfc.nb() != block_band)
        throw std::runtime_error(
            "k-BLACS eigenvector target must use the permanent capped rectangular layout");
}

librpa_int::ComplexMatrix &prepare_local_wfc(
    librpa_int::MeanField &mf, const librpa_int::ArrayDesc &desc_wfc,
    const int ispin, const int ispinor, const int ik)
{
    auto &wfc = mf.get_eigenvectors()[ispin][ispinor][ik];
    if (wfc.nr != desc_wfc.n_loc() || wfc.nc != desc_wfc.m_loc())
        wfc.create(desc_wfc.n_loc(), desc_wfc.m_loc(), false);
    return wfc;
}

int target_from_source(const int ik_source, const int n_source_kpoints,
                       const int n_target_kpoints,
                       const std::vector<int> *source_to_target_ik)
{
    if (ik_source < 0 || ik_source >= n_source_kpoints) return -2;
    if (source_to_target_ik == nullptr)
        return ik_source < n_target_kpoints ? ik_source : -2;
    return source_to_target_ik->at(static_cast<std::size_t>(ik_source));
}

int read_binary_v1_file_kblacs_2d(
    const std::string &file_path, librpa_int::MeanField &mf, const WfcShape &shape,
    const librpa_int::KPointBlacsParallelContext &kblacs_ctxt,
    const librpa_int::ArrayDesc &desc_wfc,
    const std::vector<int> *source_to_target_ik, std::vector<int> &target_hits_local)
{
    MPI_File file = MPI_FILE_NULL;
    if (MPI_File_open(kblacs_ctxt.comm_blacs_h.comm, const_cast<char *>(file_path.c_str()),
                      MPI_MODE_RDONLY, MPI_INFO_NULL, &file) != MPI_SUCCESS)
        return 1;

    std::array<std::int32_t, 6> header{};
    int io_error = 0;
    if (kblacs_ctxt.comm_blacs_h.is_root())
    {
        MPI_Status status;
        if (MPI_File_read_at(file, 0, header.data(), static_cast<int>(sizeof(header)), MPI_BYTE,
                             &status) != MPI_SUCCESS)
            io_error = 1;
    }
    MPI_Bcast(&io_error, 1, MPI_INT, 0, kblacs_ctxt.comm_blacs_h.comm);
    MPI_Bcast(header.data(), static_cast<int>(sizeof(header)), MPI_BYTE, 0,
              kblacs_ctxt.comm_blacs_h.comm);
    if (io_error)
    {
        MPI_File_close(&file);
        return 1;
    }

    const auto marker = header[0];
    const auto kind_raw = header[1];
    const int nkpoints_file = header[2];
    const int nspins_file = header[3];
    const int nstates_file = header[4];
    const int nbasis_wfc_file = header[5];
    if (marker != READER_EIGENVEC_V1_MARKER ||
        kind_raw != EIGENVEC_V1_KIND_COMPLEX_DOUBLE || nkpoints_file < 0 ||
        nspins_file != shape.nspin || nstates_file != shape.nband ||
        nbasis_wfc_file != shape.nao * shape.nsoc)
    {
        MPI_File_close(&file);
        return 1;
    }

    std::vector<KpointBlock> blocks(static_cast<std::size_t>(nkpoints_file));
    if (kblacs_ctxt.comm_blacs_h.is_root())
    {
        MPI_Offset offset = static_cast<MPI_Offset>(sizeof(header));
        MPI_Status status;
        for (auto &block : blocks)
        {
            if (MPI_File_read_at(file, offset, &block.ik, sizeof(block.ik), MPI_BYTE, &status) !=
                    MPI_SUCCESS ||
                MPI_File_read_at(file, offset + sizeof(block.ik), &block.offset,
                                 sizeof(block.offset), MPI_BYTE, &status) != MPI_SUCCESS)
            {
                io_error = 1;
                break;
            }
            offset += sizeof(block.ik) + sizeof(block.offset);
        }
    }
    MPI_Bcast(&io_error, 1, MPI_INT, 0, kblacs_ctxt.comm_blacs_h.comm);
    if (!blocks.empty())
        MPI_Bcast(blocks.data(), static_cast<int>(blocks.size() * sizeof(KpointBlock)), MPI_BYTE, 0,
                  kblacs_ctxt.comm_blacs_h.comm);
    if (io_error)
    {
        MPI_File_close(&file);
        return 1;
    }

    std::vector<std::complex<double>> dummy(1, {0.0, 0.0});
    const int local_count = desc_wfc.m_loc() * desc_wfc.n_loc();
    const int global_sizes[2]{shape.nband, shape.nao};
    const int local_sizes[2]{desc_wfc.n_loc(), desc_wfc.m_loc()};
    const int starts[2]{desc_wfc.n_loc() == 0 ? 0 : desc_wfc.indx_l2g_c(0),
                        desc_wfc.m_loc() == 0 ? 0 : desc_wfc.indx_l2g_r(0)};
    const MPI_Offset component_bytes =
        static_cast<MPI_Offset>(shape.nband) * shape.nao * sizeof(std::complex<double>);

    for (const auto &block : blocks)
    {
        const int ik_source = block.ik - 1;
        const int ik_target = target_from_source(ik_source, shape.nkpoints,
                                                 mf.get_n_kpoints(), source_to_target_ik);
        if (ik_target == -2 || ik_target >= mf.get_n_kpoints())
        {
            MPI_File_close(&file);
            return 1;
        }
        if (ik_target < 0 || !kblacs_ctxt.owns_kpoint(ik_target)) continue;
        if (kblacs_ctxt.comm_blacs_h.is_root())
            ++target_hits_local.at(static_cast<std::size_t>(ik_target));

        for (int ispin = 0; ispin != shape.nspin; ++ispin)
        {
            for (int ispinor = 0; ispinor != shape.nsoc; ++ispinor)
            {
                auto &wfc = prepare_local_wfc(mf, desc_wfc, ispin, ispinor, ik_target);
                MPI_Datatype filetype = MPI_C_DOUBLE_COMPLEX;
                bool free_filetype = false;
                if (local_count > 0)
                {
                    MPI_Type_create_subarray(2, global_sizes, local_sizes, starts, MPI_ORDER_C,
                                             MPI_C_DOUBLE_COMPLEX, &filetype);
                    MPI_Type_commit(&filetype);
                    free_filetype = true;
                }
                const MPI_Offset component =
                    static_cast<MPI_Offset>(ispin * shape.nsoc + ispinor);
                const MPI_Offset displacement =
                    static_cast<MPI_Offset>(block.offset) + component * component_bytes;
                MPI_File_set_view(file, displacement, MPI_C_DOUBLE_COMPLEX, filetype,
                                  const_cast<char *>("native"), MPI_INFO_NULL);
                MPI_Status status;
                auto *dst = wfc.c == nullptr ? dummy.data() : wfc.c;
                if (MPI_File_read_all(file, dst, local_count, MPI_C_DOUBLE_COMPLEX, &status) !=
                    MPI_SUCCESS)
                    io_error = 1;
                if (free_filetype) MPI_Type_free(&filetype);
            }
        }
    }
    MPI_File_close(&file);
    MPI_Allreduce(MPI_IN_PLACE, &io_error, 1, MPI_INT, MPI_MAX,
                  kblacs_ctxt.comm_blacs_h.comm);
    return io_error;
}

int read_legacy_text_file_kblacs_2d(
    const std::string &file_path, librpa_int::MeanField &mf, const WfcShape &shape,
    const librpa_int::KPointBlacsParallelContext &kblacs_ctxt,
    const librpa_int::ArrayDesc &desc_wfc,
    const std::vector<int> *source_to_target_ik, std::vector<int> &target_hits_local,
    const LegacyTextWfcOrder text_order)
{
    std::ifstream infile(file_path);
    if (!infile.good()) return 1;

    std::string rvalue, ivalue, kstr;
    while (infile >> kstr)
    {
        const int ik_source = std::stoi(kstr) - 1;
        const int ik_target = target_from_source(ik_source, shape.nkpoints,
                                                 mf.get_n_kpoints(), source_to_target_ik);
        if (ik_target == -2 || ik_target >= mf.get_n_kpoints()) return 1;
        const bool keep = ik_target >= 0 && kblacs_ctxt.owns_kpoint(ik_target);
        if (keep && kblacs_ctxt.comm_blacs_h.is_root())
            ++target_hits_local.at(static_cast<std::size_t>(ik_target));

        const auto store = [&](const int ispin, const int ispinor, const int ib, const int iw)
        {
            if (!keep) return;
            const int iloc = desc_wfc.indx_g2l_r(iw);
            const int jloc = desc_wfc.indx_g2l_c(ib);
            if (iloc < 0 || jloc < 0) return;
            auto &wfc = prepare_local_wfc(mf, desc_wfc, ispin, ispinor, ik_target);
            wfc(jloc, iloc) = {std::stod(rvalue), std::stod(ivalue)};
        };

        if (text_order == LegacyTextWfcOrder::BasisSpinorBandSpin)
        {
            for (int iw = 0; iw != shape.nao; ++iw)
                for (int ispinor = 0; ispinor != shape.nsoc; ++ispinor)
                    for (int ib = 0; ib != shape.nband; ++ib)
                        for (int ispin = 0; ispin != shape.nspin; ++ispin)
                        {
                            if (!(infile >> rvalue >> ivalue)) return 1;
                            store(ispin, ispinor, ib, iw);
                        }
        }
        else
        {
            for (int ispin = 0; ispin != shape.nspin; ++ispin)
                for (int iwfc = 0; iwfc != shape.nao * shape.nsoc; ++iwfc)
                {
                    const int ispinor = iwfc % shape.nsoc;
                    const int iw = iwfc / shape.nsoc;
                    for (int ib = 0; ib != shape.nband; ++ib)
                    {
                        if (!(infile >> rvalue >> ivalue)) return 1;
                        store(ispin, ispinor, ib, iw);
                    }
                }
        }
    }
    return 0;
}

}  // namespace

int read_eigenvector(ReaderContext &ctx, const std::string &dir_path)
{
    const WfcShape shape{ctx.state.n_spins,   ctx.state.n_spinor,
                         ctx.state.n_states,  ctx.state.n_basis_ao,
                         ctx.state.n_kpoints, ctx.params.use_spinor_wfc};

    int files_read = 0;
    int version_first = -1;
    for (const auto &file_path : eigenvector_files(ctx, dir_path))
    {
        librpa_int::require_readable_file(file_path);
        const int version = check_KS_file_version(file_path);
        if (version_first < 0)
        {
            version_first = version;
            librpa_int::global::lib_printf_root("KS eigenvector reader: %s\n",
                                                version == 1 ? "binary v1" : "legacy text");
        }
        else
        {
            if (version != version_first)
                throw std::runtime_error(
                    "Versions across KS eigenvector files are inconsistent! "
                    "Version of first file " + std::to_string(version_first) +
                    " vs " + std::to_string(version) + " of " + file_path);
        }
        int ret = 1;
        switch (version)
        {
            case 0:
                ret = read_legacy_text_file(
                    file_path, shape, [&](int ik) { return selected_input_ik(ctx, ik); },
                    [&](const int ik, const std::vector<double> &re, const std::vector<double> &im)
                    { set_wfc(ctx, shape, ik, re, im); });
                break;
            case 1:
                ret = read_binary_v1_file(
                    file_path, shape, [&](int ik) { return selected_input_ik(ctx, ik); },
                    [&](const int ik, const std::vector<std::complex<double>> &block_data)
                    { set_wfc_packed(ctx, shape, ik, block_data); });
                break;
        }
        if (ret != 0) return ret;
        ++files_read;
    }
    return files_read == 0 ? -1 : 0;
}

int read_eigenvector_kblacs_2d(ReaderContext &ctx,
    const std::string &dir_path, librpa_int::MeanField &mf, const bool use_spinor_wfc,
    const librpa_int::KPointBlacsParallelContext &kblacs_ctxt,
    const librpa_int::ArrayDesc &desc_wfc,
    const std::vector<int> *source_to_target_ik, const LegacyTextWfcOrder text_order)
{
    validate_kblacs_wfc_layout(mf, kblacs_ctxt, desc_wfc);
    const auto desc_wfc_io =
        kblacs_ctxt.create_array_desc(mf.get_n_aos(), mf.get_n_states());
    if (!desc_wfc_io.is_row_consec() || !desc_wfc_io.is_col_consec())
        throw std::runtime_error(
            "temporary direct eigenvector input descriptor must be contiguous");
    const int n_source_kpoints =
        source_to_target_ik == nullptr ? mf.get_n_kpoints()
                                       : static_cast<int>(source_to_target_ik->size());
    if (n_source_kpoints <= 0) return 1;

    std::vector<int> expected_hits(static_cast<std::size_t>(mf.get_n_kpoints()), 0);
    for (int ik_source = 0; ik_source != n_source_kpoints; ++ik_source)
    {
        const int ik_target = target_from_source(ik_source, n_source_kpoints,
                                                 mf.get_n_kpoints(), source_to_target_ik);
        if (ik_target < -1 || ik_target >= mf.get_n_kpoints()) return 1;
        if (ik_target >= 0) ++expected_hits[static_cast<std::size_t>(ik_target)];
    }
    for (const int hits : expected_hits)
        if (hits != 1) return 1;

    const WfcShape shape{mf.get_n_spins(), mf.get_n_spinor(), mf.get_n_states(),
                         mf.get_n_aos(), n_source_kpoints, use_spinor_wfc};
    mf.get_eigenvectors().clear();
    std::vector<int> target_hits_local(static_cast<std::size_t>(mf.get_n_kpoints()), 0);
    int files_read = 0;
    int version_first = -1;
    for (const auto &file_path : eigenvector_files(ctx, dir_path))
    {
        librpa_int::require_readable_file(file_path);
        const int version = check_KS_file_version(file_path);
        if (version_first < 0)
        {
            version_first = version;
            librpa_int::global::lib_printf_root(
                "KS eigenvector reader: %s, direct k-BLACS 2D input\n",
                version == 1 ? "binary v1 MPI-IO" : "legacy text streaming");
        }
        else if (version != version_first)
        {
            throw std::runtime_error("Versions across KS eigenvector files are inconsistent");
        }

        const int ret = version == 1
                            ? read_binary_v1_file_kblacs_2d(
                                  file_path, mf, shape, kblacs_ctxt, desc_wfc_io,
                                  source_to_target_ik, target_hits_local)
                            : read_legacy_text_file_kblacs_2d(
                                  file_path, mf, shape, kblacs_ctxt, desc_wfc_io,
                                  source_to_target_ik, target_hits_local, text_order);
        if (ret != 0) return ret;
        ++files_read;
    }
    if (files_read == 0) return -1;

    std::vector<int> target_hits(target_hits_local.size(), 0);
    MPI_Allreduce(target_hits_local.data(), target_hits.data(),
                  static_cast<int>(target_hits.size()), MPI_INT, MPI_SUM,
                  kblacs_ctxt.comm_global_h.comm);
    for (std::size_t ik = 0; ik != target_hits.size(); ++ik)
        if (target_hits[ik] != expected_hits[ik]) return 1;
    librpa_int::redistribute_meanfield_eigvecs_kblacs(
        mf, kblacs_ctxt, desc_wfc_io, desc_wfc, "direct SCF");
    return 0;
}

int read_eigenvector(ReaderContext &ctx, const std::string &dir_path, librpa_int::MeanField &mf, bool use_spinor_wfc,
                     const std::vector<int> *iks_selected)
{
    const WfcShape shape{mf.get_n_spins(), mf.get_n_spinor(),  mf.get_n_states(),
                         mf.get_n_aos(),   mf.get_n_kpoints(), use_spinor_wfc};

    int files_read = 0;
    bool printed_reader_version = false;
    for (const auto &file_path : eigenvector_files(ctx, dir_path))
    {
        librpa_int::require_readable_file(file_path);
        const int version = check_KS_file_version(file_path);
        if (!printed_reader_version)
        {
            librpa_int::global::lib_printf_root("KS eigenvector reader: %s\n",
                                                version == 1 ? "binary v1" : "legacy text");
            printed_reader_version = true;
        }
        const auto should_read_ik = [&](const int ik) { return selected_ik(iks_selected, ik); };
        const auto store_ik =
            [&](const int ik, const std::vector<double> &re, const std::vector<double> &im)
        { set_meanfield_wfc(mf, shape, ik, re, im); };
        const auto ret =
            version == 1
                ? read_binary_v1_file(
                      file_path, shape, should_read_ik,
                      [&](const int ik, const std::vector<std::complex<double>> &block_data)
                      { set_meanfield_wfc_packed(mf, shape, ik, block_data); })
                : read_legacy_text_file(file_path, shape, should_read_ik, store_ik);
        if (ret != 0) return ret;
        ++files_read;
    }
    return files_read == 0 ? -1 : 0;
}

int read_eigenvector(ReaderContext &ctx, const std::string &dir_path, librpa_int::MeanField &mf, bool use_spinor_wfc,
                     const std::vector<int> &source_to_target_ik,
                     const std::vector<int> *source_iks_selected,
                     const LegacyTextWfcOrder text_order)
{
    const int n_target_kpoints = mf.get_n_kpoints();
    const int n_source_kpoints = static_cast<int>(source_to_target_ik.size());
    if (n_source_kpoints <= 0)
    {
        return 1;
    }
    for (const int ik_target : source_to_target_ik)
    {
        if (ik_target < -1 || ik_target >= n_target_kpoints)
        {
            return 1;
        }
    }

    const WfcShape shape{mf.get_n_spins(), mf.get_n_spinor(),  mf.get_n_states(),
                         mf.get_n_aos(),   n_source_kpoints,   use_spinor_wfc};

    std::vector<char> target_hits(static_cast<std::size_t>(n_target_kpoints), 0);
    const auto source_selected = [&](const int ik_source)
    {
        if (ik_source < 0 || ik_source >= n_source_kpoints) return false;
        if (source_to_target_ik[ik_source] < 0) return false;
        return selected_ik(source_iks_selected, ik_source);
    };
    const auto target_from_source = [&](const int ik_source)
    {
        return source_to_target_ik.at(static_cast<std::size_t>(ik_source));
    };

    int files_read = 0;
    bool printed_reader_version = false;
    for (const auto &file_path : eigenvector_files(ctx, dir_path))
    {
        const int version = check_KS_file_version(file_path);
        if (!printed_reader_version)
        {
            librpa_int::global::lib_printf_root("KS eigenvector reader: %s\n",
                                                version == 1 ? "binary v1" : "legacy text");
            printed_reader_version = true;
        }

        const auto store_text =
            [&](const int ik_source, const std::vector<double> &re,
                const std::vector<double> &im)
        {
            const int ik_target = target_from_source(ik_source);
            set_meanfield_wfc(mf, shape, ik_target, re, im);
            target_hits.at(static_cast<std::size_t>(ik_target)) = 1;
        };
        const auto store_packed =
            [&](const int ik_source, const std::vector<std::complex<double>> &block_data)
        {
            const int ik_target = target_from_source(ik_source);
            set_meanfield_wfc_packed(mf, shape, ik_target, block_data);
            target_hits.at(static_cast<std::size_t>(ik_target)) = 1;
        };

        const auto ret =
            version == 1
                ? read_binary_v1_file(file_path, shape, source_selected, store_packed)
                : read_legacy_text_file(file_path, shape, source_selected, store_text,
                                        text_order);
        if (ret != 0) return ret;
        ++files_read;
    }

    if (files_read == 0) return -1;
    if (source_iks_selected == nullptr)
    {
        for (const auto hit : target_hits)
        {
            if (!hit) return 1;
        }
    }
    return 0;
}

} // namespace librpa::reader
