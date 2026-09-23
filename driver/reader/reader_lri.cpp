#include "reader_lri.h"

#include <dirent.h>

#include <algorithm>
#include <cassert>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "reader_context.h"
#include "../../src/io/fs.h"
#include "../../src/io/global_io.h"
#include "../../src/io/stl_io_helper.h"
#include "../../src/mpi/global_mpi.h"
#include "../../src/utils/error.h"
#include "../../src/utils/profiler.h"

namespace librpa::reader
{

#define READER_LRICOEF_V1_MARKER -10267453

using std::ifstream;
using std::string;
using librpa_int::atpair_t;
using librpa_int::matrix;

static bool starts_with_prefix(const std::string &text, const std::string &prefix)
{
    return text.rfind(prefix, 0) == 0;
}

static std::string basename_from_path(const std::string &path)
{
    const auto pos = path.find_last_of("/\\");
    if (pos == std::string::npos)
    {
        return path;
    }
    return path.substr(pos + 1);
}

static void assert_distinct_Cs_prefixes(ReaderContext &ctx)
{
    const auto &full_prefix = ctx.params.prefix_lri_coeff;
    const auto &shrink_prefix = ctx.params.prefix_lri_coeff_shrink;
    if (!full_prefix.empty() && full_prefix == shrink_prefix)
    {
        throw LIBRPA_RUNTIME_ERROR(
            "prefix_lri_coeff and prefix_lri_coeff_shrink must be distinct; "
            "reader-v1 full Cs and shrink Cs are separate file families");
    }
}

static std::vector<std::string> discover_Cs_files_for_keyword(ReaderContext &ctx, const std::string &dir_path,
                                                             const std::string &keyword)
{
    assert_distinct_Cs_prefixes(ctx);

    auto files = librpa_int::discover_files_with_prefix(dir_path, keyword);
    const auto &full_prefix = ctx.params.prefix_lri_coeff;
    const auto &shrink_prefix = ctx.params.prefix_lri_coeff_shrink;
    std::string excluded_prefix;
    if (keyword == full_prefix && starts_with_prefix(shrink_prefix, full_prefix))
    {
        excluded_prefix = shrink_prefix;
    }
    else if (keyword == shrink_prefix && starts_with_prefix(full_prefix, shrink_prefix))
    {
        excluded_prefix = full_prefix;
    }

    if (!excluded_prefix.empty())
    {
        files.erase(std::remove_if(files.begin(), files.end(),
                                   [&excluded_prefix](const std::string &path)
                                   {
                                       return starts_with_prefix(basename_from_path(path),
                                                                 excluded_prefix);
                                   }),
                    files.end());
    }
    return files;
}

static bool target_shrink_aux_basis_for_keyword(ReaderContext &ctx, const std::string &keyword)
{
    if (keyword == ctx.params.prefix_lri_coeff_shrink)
        return true;
    return false;
}

static const std::vector<size_t> &target_aux_basis_sizes_for_keyword(ReaderContext &ctx, const std::string &keyword)
{
    if (target_shrink_aux_basis_for_keyword(ctx, keyword))
        return ctx.state.nbs_aux_shrink;
    return ctx.state.nbs_aux;
}

static void set_lri_coeff(ReaderContext &ctx, const std::string &keyword, LibrpaParallelRouting routing,
                                 int I, int J, int nbasis_i, int nbasis_j, int naux_mu,
                                 const int R[3], const double *Cs_in)
{
    const int shrink_aux = keyword == ctx.params.prefix_lri_coeff_shrink ? 1 : 0;
    ctx.h.set_lri_coeff(routing, I, J, nbasis_i, nbasis_j, naux_mu, R, Cs_in, shrink_aux);
}

static bool get_Cs_binary_data_size(int n_i, int n_j, int n_mu, std::streamoff &data_size)
{
    if (n_i <= 0 || n_j <= 0 || n_mu <= 0) return false;

    const auto max_count = static_cast<unsigned long long>(
        std::numeric_limits<std::streamoff>::max() / sizeof(double));
    auto count = static_cast<unsigned long long>(n_i);
    const auto n_j_ull = static_cast<unsigned long long>(n_j);
    const auto n_mu_ull = static_cast<unsigned long long>(n_mu);

    if (count > max_count / n_j_ull) return false;
    count *= n_j_ull;
    if (count > max_count / n_mu_ull) return false;
    count *= n_mu_ull;

    data_size = static_cast<std::streamoff>(count * sizeof(double));
    return true;
}

static_assert(sizeof(int) == sizeof(std::int32_t),
              "Cs binary readers expect native int to be int32");

namespace
{

constexpr std::streamoff CS_BINARY_LEGACY_HEADER_SIZE =
    3 * static_cast<std::streamoff>(sizeof(int));
constexpr std::streamoff CS_BINARY_V1_HEADER_BASE_SIZE =
    3 * static_cast<std::streamoff>(sizeof(int)) +
    2 * static_cast<std::streamoff>(sizeof(std::int64_t));
constexpr std::streamoff CS_BINARY_LEGACY_BLOCK_HEADER_SIZE =
    8 * static_cast<std::streamoff>(sizeof(int));
constexpr std::size_t CS_BINARY_V1_BLOCK_RECORD_SIZE =
    5 * sizeof(std::int32_t) + sizeof(double) + sizeof(std::int64_t);
constexpr std::size_t CS_BINARY_V1_MAX_COLLECTED_READ_BYTES =
    256ULL * 1024ULL * 1024ULL;

struct CsBinaryV1Record
{
    int ia1 = 0;
    int ia2 = 0;
    int R[3] = {0, 0, 0};
    double max_abs = 0.0;
    std::int64_t offset = 0;
};

struct CsBinaryV1ReadTask
{
    std::string file_path;
    std::size_t block_id = 0;
    CsBinaryV1Record block;
    std::streamoff nbytes = 0;
    int n_i = 0;
    int n_j = 0;
    int n_mu = 0;
};

struct CsBinaryV1CollectedRead
{
    std::string file_path;
    std::streamoff offset = 0;
    std::streamoff nbytes = 0;
    std::vector<CsBinaryV1ReadTask> tasks;
};

bool get_Cs_binary_v1_header_size(const std::int64_t nblocks_max,
                                  std::streamoff &header_size)
{
    if (nblocks_max < 0 ||
        static_cast<unsigned long long>(nblocks_max) >
            static_cast<unsigned long long>(std::numeric_limits<std::streamoff>::max()))
    {
        return false;
    }
    const auto record_size = static_cast<std::streamoff>(CS_BINARY_V1_BLOCK_RECORD_SIZE);
    const auto max_header_size = std::numeric_limits<std::streamoff>::max();
    if (static_cast<std::streamoff>(nblocks_max) >
        (max_header_size - CS_BINARY_V1_HEADER_BASE_SIZE) / record_size)
        return false;
    header_size =
        CS_BINARY_V1_HEADER_BASE_SIZE + static_cast<std::streamoff>(nblocks_max) * record_size;
    return true;
}

bool parse_Cs_binary_v1_header(const string &file_path,
                               std::vector<CsBinaryV1Record> &records,
                               int &natom,
                               int &ncell,
                               std::streamoff &file_size,
                               std::string *error)
{
    auto fail = [error](const std::string &message)
    {
        if (error != nullptr) *error = message;
        return false;
    };

    records.clear();
    natom = 0;
    ncell = 0;
    file_size = 0;

    librpa_int::require_readable_file(file_path);
    ifstream infile(file_path, std::ios::in | std::ios::binary);
    if (!infile.good())
    {
        return fail("Failed to open " + file_path);
    }

    infile.seekg(0, std::ios::end);
    const auto end_pos = infile.tellg();
    if (end_pos == std::streampos(-1))
    {
        return fail("Failed to determine size of binary Cs file: " + file_path);
    }
    file_size = static_cast<std::streamoff>(end_pos);
    if (file_size < CS_BINARY_V1_HEADER_BASE_SIZE)
    {
        return fail("Binary Cs file is too small: " + file_path);
    }

    int marker = 0;
    std::int64_t n_apcell_file = 0;
    std::int64_t n_apcell_file_max = 0;
    infile.seekg(0, std::ios::beg);
    infile.read((char *) &marker, sizeof(int));
    infile.read((char *) &natom, sizeof(int));
    infile.read((char *) &ncell, sizeof(int));
    infile.read((char *) &n_apcell_file, sizeof(std::int64_t));
    infile.read((char *) &n_apcell_file_max, sizeof(std::int64_t));
    if (!infile.good() || marker != READER_LRICOEF_V1_MARKER)
    {
        return fail("Invalid binary Cs v1 marker in: " + file_path);
    }
    if (!infile.good() || natom <= 0 || ncell < 0 || n_apcell_file < 0 ||
        n_apcell_file_max < n_apcell_file)
    {
        return fail("Invalid binary Cs v1 header in: " + file_path);
    }

    std::streamoff header_size = 0;
    if (!get_Cs_binary_v1_header_size(n_apcell_file_max, header_size) ||
        header_size > file_size)
    {
        return fail("Invalid binary Cs v1 header size in: " + file_path);
    }
    records.reserve(static_cast<std::size_t>(n_apcell_file));
    for (std::int64_t i = 0; i != n_apcell_file_max; ++i)
    {
        char record[CS_BINARY_V1_BLOCK_RECORD_SIZE];
        infile.read(record, CS_BINARY_V1_BLOCK_RECORD_SIZE);
        if (!infile.good())
        {
            return fail("Truncated binary Cs v1 block table in: " + file_path);
        }

        std::int32_t ints[5];
        for (int j = 0; j != 5; ++j)
        {
            std::memcpy(&ints[j], record + j * sizeof(std::int32_t), sizeof(std::int32_t));
        }

        CsBinaryV1Record block;
        block.ia1 = static_cast<int>(ints[0]);
        block.ia2 = static_cast<int>(ints[1]);
        block.R[0] = static_cast<int>(ints[2]);
        block.R[1] = static_cast<int>(ints[3]);
        block.R[2] = static_cast<int>(ints[4]);
        std::memcpy(&block.max_abs, record + 5 * sizeof(std::int32_t), sizeof(double));
        std::memcpy(&block.offset,
                    record + 5 * sizeof(std::int32_t) + sizeof(double),
                    sizeof(std::int64_t));

        if (i >= n_apcell_file)
        {
            const bool is_zero_padding =
                block.ia1 == 0 && block.ia2 == 0 &&
                block.R[0] == 0 && block.R[1] == 0 && block.R[2] == 0 &&
                block.max_abs == 0.0 && block.offset == 0;
            if (!is_zero_padding)
            {
                return fail("Nonzero padding record in binary Cs v1 block table: " +
                            file_path);
            }
            continue;
        }

        if (block.ia1 <= 0 || block.ia1 > natom ||
            block.ia2 <= 0 || block.ia2 > natom)
        {
            return fail("Invalid atom index in binary Cs v1 block table: " + file_path);
        }
        if (!std::isfinite(block.max_abs) || block.max_abs < 0.0)
        {
            return fail("Invalid max-abs value in binary Cs v1 block table: " + file_path);
        }
        if (block.offset < header_size || block.offset >= file_size)
        {
            return fail("Invalid byte offset in binary Cs v1 block table: " + file_path);
        }

        records.push_back(block);
    }
    return true;
}

bool has_Cs_binary_v1_layout(const string &file_path)
{
    std::vector<CsBinaryV1Record> records;
    int natom = 0;
    int ncell = 0;
    std::streamoff file_size = 0;
    return parse_Cs_binary_v1_header(file_path, records, natom, ncell, file_size, nullptr);
}

std::vector<CsBinaryV1Record> read_Cs_binary_v1_header_or_throw(
    const string &file_path,
    int &natom,
    int &ncell,
    std::streamoff &file_size)
{
    std::vector<CsBinaryV1Record> records;
    std::string error;
    if (!parse_Cs_binary_v1_header(file_path, records, natom, ncell, file_size, &error))
    {
        throw std::runtime_error(error);
    }
    return records;
}

void throw_unknown_Cs_reader_version(const int reader_version)
{
    throw std::logic_error("Unknown LRI coefficient reader version " +
                           std::to_string(reader_version));
}

void throw_Cs_v1_requires_reader_version(const string &file_path)
{
    throw std::logic_error(
        "LRI coefficient reader v1 file found while version_lri_reader = 0; "
        "set version_lri_reader = 1 to read: " + file_path);
}

int checked_Cs_basis_size(const std::size_t value,
                          const string &file_path,
                          const char *basis_name,
                          const int atom)
{
    if (value == 0 ||
        value > static_cast<std::size_t>(std::numeric_limits<int>::max()))
    {
        std::ostringstream ss;
        ss << file_path << ": invalid " << basis_name << " size " << value
           << " for atom " << atom;
        throw std::runtime_error(ss.str());
    }
    return static_cast<int>(value);
}

void get_Cs_binary_v1_block_dimensions(const CsBinaryV1Record &block,
                                       const std::vector<size_t> &basis_wfc,
                                       const std::vector<size_t> &basis_aux,
                                       const string &file_path,
                                       int &n_i,
                                       int &n_j,
                                       int &n_mu)
{
    const int ia1 = block.ia1 - 1;
    const int ia2 = block.ia2 - 1;
    n_i = checked_Cs_basis_size(basis_wfc[ia1], file_path, "wave-function basis", ia1);
    n_j = checked_Cs_basis_size(basis_wfc[ia2], file_path, "wave-function basis", ia2);
    n_mu = checked_Cs_basis_size(basis_aux[ia1], file_path, "auxiliary basis", ia1);
}

void validate_Cs_binary_v1_blocks(const string &file_path,
                                  const int natom,
                                  const std::streamoff file_size,
                                  const std::vector<CsBinaryV1Record> &records,
                                  const std::vector<size_t> &basis_wfc,
                                  const std::vector<size_t> &basis_aux)
{
    if (basis_wfc.empty())
    {
        throw std::runtime_error(file_path +
                                 ": Cs reader v1 requires wave-function basis information");
    }
    if (basis_aux.empty())
    {
        throw std::runtime_error(file_path +
                                 ": Cs reader v1 requires auxiliary basis information");
    }
    if (basis_wfc.size() != static_cast<std::size_t>(natom) ||
        basis_aux.size() != static_cast<std::size_t>(natom))
    {
        std::ostringstream ss;
        ss << file_path << ": Cs reader v1 atom count " << natom
           << " does not match loaded basis atom counts wfc=" << basis_wfc.size()
           << ", aux=" << basis_aux.size();
        throw std::runtime_error(ss.str());
    }

    std::vector<std::pair<std::streamoff, std::streamoff>> ranges;
    ranges.reserve(records.size());
    for (const auto &block: records)
    {
        int n_i = 0;
        int n_j = 0;
        int n_mu = 0;
        get_Cs_binary_v1_block_dimensions(
            block, basis_wfc, basis_aux, file_path, n_i, n_j, n_mu);

        std::streamoff data_size = 0;
        if (!get_Cs_binary_data_size(n_i, n_j, n_mu, data_size) ||
            block.offset < 0 ||
            block.offset > std::numeric_limits<std::streamoff>::max() ||
            data_size > file_size - static_cast<std::streamoff>(block.offset))
        {
            std::ostringstream ss;
            ss << file_path << ": invalid Cs reader v1 payload offset "
               << block.offset << " for atom pair (" << block.ia1 << ", "
               << block.ia2 << ")";
            throw std::runtime_error(ss.str());
        }
        const auto begin = static_cast<std::streamoff>(block.offset);
        ranges.emplace_back(begin, begin + data_size);
    }

    std::sort(ranges.begin(), ranges.end());
    for (std::size_t i = 1; i != ranges.size(); ++i)
    {
        if (ranges[i].first < ranges[i - 1].second)
        {
            throw std::runtime_error(file_path + ": overlapping Cs reader v1 payload blocks");
        }
    }
}

std::vector<CsBinaryV1ReadTask> make_Cs_binary_v1_read_tasks(
    const std::vector<string> &files,
    const double threshold,
    const std::vector<size_t> &basis_wfc,
    const std::vector<size_t> &basis_aux)
{
    std::vector<CsBinaryV1ReadTask> tasks;
    for (const auto &file_path: files)
    {
        int natom = 0;
        int ncell = 0;
        std::streamoff file_size = 0;
        const auto records = read_Cs_binary_v1_header_or_throw(
            file_path, natom, ncell, file_size);
        validate_Cs_binary_v1_blocks(
            file_path, natom, file_size, records, basis_wfc, basis_aux);

        std::vector<CsBinaryV1ReadTask> file_tasks;
        file_tasks.reserve(records.size());
        for (std::size_t i = 0; i != records.size(); ++i)
        {
            const auto &block = records[i];
            if (block.max_abs < threshold) continue;

            int n_i = 0;
            int n_j = 0;
            int n_mu = 0;
            get_Cs_binary_v1_block_dimensions(
                block, basis_wfc, basis_aux, file_path, n_i, n_j, n_mu);

            std::streamoff data_size = 0;
            if (!get_Cs_binary_data_size(n_i, n_j, n_mu, data_size))
            {
                throw std::runtime_error("Invalid Cs reader v1 block dimensions in: " +
                                         file_path);
            }
            file_tasks.push_back({file_path, i, block, data_size, n_i, n_j, n_mu});
        }

        std::sort(file_tasks.begin(), file_tasks.end(),
                  [](const auto &a, const auto &b)
                  {
                      if (a.block.offset != b.block.offset)
                          return a.block.offset < b.block.offset;
                      return a.block_id < b.block_id;
                  });
        tasks.insert(tasks.end(), file_tasks.begin(), file_tasks.end());
    }
    return tasks;
}

std::uint64_t checked_task_nbytes(const std::streamoff nbytes,
                                  const std::string &file_path)
{
    if (nbytes < 0 ||
        static_cast<unsigned long long>(nbytes) >
            static_cast<unsigned long long>(std::numeric_limits<std::uint64_t>::max()))
    {
        throw std::runtime_error(file_path + ": Cs reader v1 task byte size is out of range");
    }
    return static_cast<std::uint64_t>(nbytes);
}

int owner_rank_from_cumulative_begin(const std::uint64_t cumulative_begin,
                                     const std::uint64_t total_bytes,
                                     const int nprocs)
{
    if (total_bytes == 0) return 0;
    const long double scaled =
        static_cast<long double>(cumulative_begin) * nprocs /
        static_cast<long double>(total_bytes);
    int owner = static_cast<int>(scaled);
    if (owner < 0) owner = 0;
    if (owner >= nprocs) owner = nprocs - 1;
    return owner;
}

std::vector<CsBinaryV1ReadTask> select_Cs_binary_v1_tasks_for_rank(
    const std::vector<CsBinaryV1ReadTask> &tasks,
    const int myid,
    const int nprocs)
{
    std::uint64_t total_bytes = 0;
    for (const auto &task: tasks)
    {
        const auto nbytes = checked_task_nbytes(task.nbytes, task.file_path);
        if (total_bytes > std::numeric_limits<std::uint64_t>::max() - nbytes)
        {
            throw std::runtime_error("Cs reader v1 total scheduled bytes are out of range");
        }
        total_bytes += nbytes;
    }

    std::vector<CsBinaryV1ReadTask> selected;
    std::uint64_t cumulative_begin = 0;
    for (const auto &task: tasks)
    {
        const auto nbytes = checked_task_nbytes(task.nbytes, task.file_path);
        if (owner_rank_from_cumulative_begin(cumulative_begin, total_bytes, nprocs) == myid)
        {
            selected.push_back(task);
        }
        cumulative_begin += nbytes;
    }
    return selected;
}

std::size_t checked_streamoff_to_size(const std::streamoff value,
                                      const std::string &file_path)
{
    if (value < 0 ||
        static_cast<unsigned long long>(value) >
            static_cast<unsigned long long>(std::numeric_limits<std::size_t>::max()))
    {
        throw std::runtime_error(file_path + ": Cs reader v1 byte count is out of range");
    }
    return static_cast<std::size_t>(value);
}

std::vector<CsBinaryV1CollectedRead> collect_Cs_binary_v1_reads(
    const std::vector<CsBinaryV1ReadTask> &tasks)
{
    std::vector<CsBinaryV1CollectedRead> collected;
    for (const auto &task: tasks)
    {
        const auto task_begin = static_cast<std::streamoff>(task.block.offset);
        const auto task_end = task_begin + task.nbytes;
        if (task_end < task_begin)
        {
            throw std::runtime_error(task.file_path + ": Cs reader v1 task range is invalid");
        }

        bool need_new_read = collected.empty() ||
            task.file_path != collected.back().file_path ||
            task_begin < collected.back().offset + collected.back().nbytes;

        if (!need_new_read)
        {
            const auto read_begin = collected.back().offset;
            const auto read_nbytes = task_end - read_begin;
            need_new_read =
                read_nbytes < 0 ||
                static_cast<unsigned long long>(read_nbytes) >
                    static_cast<unsigned long long>(CS_BINARY_V1_MAX_COLLECTED_READ_BYTES);
        }

        if (need_new_read)
        {
            collected.push_back({task.file_path, task_begin, task.nbytes, {task}});
        }
        else
        {
            auto &read = collected.back();
            read.nbytes = task_end - read.offset;
            read.tasks.push_back(task);
        }
    }
    return collected;
}

void read_Cs_binary_v1_tasks(ReaderContext &ctx, const std::vector<CsBinaryV1ReadTask> &tasks,
                             const std::string &keyword)
{
    const auto collected_reads = collect_Cs_binary_v1_reads(tasks);
    for (const auto &read: collected_reads)
    {
        ifstream infile(read.file_path, std::ios::in | std::ios::binary);
        if (!infile.good())
        {
            throw std::logic_error("Failed to open " + read.file_path);
        }

        std::vector<char> buffer(checked_streamoff_to_size(read.nbytes, read.file_path));
        infile.seekg(read.offset, std::ios::beg);
        infile.read(buffer.data(), read.nbytes);
        if (!infile.good())
        {
            throw std::runtime_error("Failed to read Cs reader v1 payload in: " +
                                     read.file_path);
        }

        for (const auto &task: read.tasks)
        {
            std::shared_ptr<matrix> cs_ptr = std::make_shared<matrix>();
            cs_ptr->create(task.n_i * task.n_j, task.n_mu);
            const auto buffer_offset =
                checked_streamoff_to_size(
                    static_cast<std::streamoff>(task.block.offset) - read.offset,
                    task.file_path);
            std::memcpy(cs_ptr->c, buffer.data() + buffer_offset,
                        checked_streamoff_to_size(task.nbytes, task.file_path));

            const int ia1 = task.block.ia1 - 1;
            const int ia2 = task.block.ia2 - 1;
            int R[3] = {task.block.R[0], task.block.R[1], task.block.R[2]};
            set_lri_coeff(ctx, keyword, ctx.opts.parallel_routing, ia1, ia2,
                                 task.n_i, task.n_j, task.n_mu, R, cs_ptr->c);
        }
    }
}

} // namespace

int detect_Cs_reader_version(ReaderContext &ctx, const string &dir_path, const string keyword)
{
    const auto files = discover_Cs_files_for_keyword(ctx, dir_path, keyword);
    if (files.empty())
    {
        throw std::logic_error("No LRI coefficient files found with prefix " + keyword);
    }

    ifstream infile(files.front(), std::ios::in | std::ios::binary);
    if (!infile.good())
    {
        throw std::logic_error("Failed to open " + files.front());
    }

    int marker_or_natom = 0;
    infile.read((char *) &marker_or_natom, sizeof(int));
    if (infile.good())
    {
        if (marker_or_natom >= 0) return 0;
        if (marker_or_natom == READER_LRICOEF_V1_MARKER) return 1;

        throw std::logic_error(
            files.front() + ": unknown LRI coefficient reader marker " +
            std::to_string(marker_or_natom));
    }

    int text_value = 0;
    ifstream text_infile(files.front(), std::ios::in);
    if (text_infile >> text_value && text_value >= 0)
    {
        return 0;
    }

    throw std::logic_error("Failed to detect LRI coefficient reader version from " +
                           files.front());
}

static bool has_Cs_binary_layout(const string &file_path)
{
    librpa_int::require_readable_file(file_path);
    ifstream infile(file_path, std::ios::in | std::ios::binary);
    if (!infile.good())
    {
        throw std::logic_error("Failed to open " + file_path);
    }

    infile.seekg(0, std::ios::end);
    const auto end_pos = infile.tellg();
    if (end_pos == std::streampos(-1)) return false;

    const auto file_size = static_cast<std::streamoff>(end_pos);
    const auto header_size = CS_BINARY_LEGACY_HEADER_SIZE;
    const auto block_header_size = CS_BINARY_LEGACY_BLOCK_HEADER_SIZE;
    if (file_size < header_size) return false;

    infile.seekg(0, std::ios::beg);
    int natom = 0;
    int ncell = 0;
    int n_apcell_file = 0;
    infile.read((char *) &natom, sizeof(int));
    infile.read((char *) &ncell, sizeof(int));
    infile.read((char *) &n_apcell_file, sizeof(int));
    if (!infile.good() || natom <= 0 || ncell < 0 || n_apcell_file < 0) return false;

    std::streamoff pos = header_size;
    if (n_apcell_file > (file_size - pos) / block_header_size) return false;

    for (int i = 0; i < n_apcell_file; i++)
    {
        if (file_size - pos < block_header_size) return false;

        int dims[8];
        infile.seekg(pos, std::ios::beg);
        infile.read((char *) &dims[0], 8 * sizeof(int));
        if (!infile.good()) return false;
        pos += block_header_size;

        if (dims[0] <= 0 || dims[0] > natom || dims[1] <= 0 || dims[1] > natom)
            return false;

        std::streamoff data_size = 0;
        if (!get_Cs_binary_data_size(dims[5], dims[6], dims[7], data_size)) return false;
        if (file_size - pos < data_size) return false;

        pos += data_size;
    }

    return pos == file_size;
}

static bool has_binary_control_bytes(const string &file_path)
{
    librpa_int::require_readable_file(file_path);
    ifstream infile(file_path, std::ios::in | std::ios::binary);
    if (!infile.good())
    {
        throw std::logic_error("Failed to open " + file_path);
    }

    char buf[256];
    infile.read(buf, sizeof(buf));
    const auto nread = infile.gcount();
    for (std::streamsize i = 0; i < nread; i++)
    {
        const auto c = static_cast<unsigned char>(buf[i]);
        if (c == '\0' || (!std::isprint(c) && !std::isspace(c))) return true;
    }
    return false;
}

static bool has_Cs_text_header(const string &file_path)
{
    librpa_int::require_readable_file(file_path);
    ifstream infile(file_path, std::ios::in);
    if (!infile.good())
    {
        throw std::logic_error("Failed to open " + file_path);
    }

    int natom = 0;
    int ncell = 0;
    if (!(infile >> natom >> ncell)) return false;
    return natom > 0 && ncell >= 0;
}

//! Check if Cs data file is in ASCII text or unformatted binary format
static bool check_Cs_file_binary(const string &file_path)
{
    // Binary Cs headers start with three native ints.  The old text-first probe could
    // mis-detect a binary file when the first byte of natom happened to be an ASCII digit.
    if (has_Cs_binary_v1_layout(file_path)) return true;
    if (has_Cs_binary_layout(file_path)) return true;
    if (has_binary_control_bytes(file_path)) return true;
    return !has_Cs_text_header(file_path);
}

static size_t handle_Cs_file(ReaderContext &ctx, const string &file_path, double threshold,
                             const std::vector<atpair_t> &local_atpair,
                             const string &keyword)
{
    using namespace std;

    librpa_int::require_readable_file(file_path);
    set<size_t> loc_atp_index;
    for (auto &lap : local_atpair)
    {
        loc_atp_index.insert(lap.first);
        loc_atp_index.insert(lap.second);
    }
    // cout<<"READING Cs from file: "<<file_path<<"  Cs_first_size:
    // "<<loc_atp_index.size()<<endl;
    // map<size_t,map<size_t,map<Vector3_Order<int>,std::shared_ptr<matrix>>>> Cs_m;
    size_t cs_discard = 0;
    string natom_s, ncell_s, ia1_s, ia2_s, ic_1, ic_2, ic_3, i_s, j_s, mu_s, Cs_ele;
    int R[3];
    ifstream infile;
    infile.open(file_path);
    infile >> natom_s >> ncell_s;
    // int natom = stoi(natom_s);
    // int ncell = stoi(ncell_s);

    /* cout<<"  Natom  Ncell  "<<natom<<"  "<<ncell<<endl; */
    // for(int loop=0;loop!=natom*natom*ncell;loop++)
    while (infile.peek() != EOF)
    {
        infile >> ia1_s >> ia2_s >> ic_1 >> ic_2 >> ic_3 >> i_s;
        if (infile.peek() == EOF) break;
        // cout << " ia1_s,ia2_s: " << ia1_s << "  " << ia2_s << endl;
        infile >> j_s >> mu_s;
        // cout<<ic_1<<mu_s<<endl;
        int ia1 = stoi(ia1_s) - 1;
        int ia2 = stoi(ia2_s) - 1;
        R[0] = stoi(ic_1);
        R[1] = stoi(ic_2);
        R[2] = stoi(ic_3);
        int n_i = stoi(i_s);
        int n_j = stoi(j_s);
        int n_mu = stoi(mu_s);

        // cout<< ia1<<ia2<<box<<endl;
        shared_ptr<matrix> cs_ptr = make_shared<matrix>();
        cs_ptr->create(n_i * n_j, n_mu);
        // cout<<cs_ptr->nr<<cs_ptr->nc<<endl;

        for (int i = 0; i != n_i; i++)
            for (int j = 0; j != n_j; j++)
                for (int mu = 0; mu != n_mu; mu++)
                {
                    infile >> Cs_ele;
                    (*cs_ptr)(i *n_j + j, mu) = stod(Cs_ele);
                    // if (i == j)
                    // {
                    //     (*cs_ptr)(i * n_j + j, mu) = 1.0;
                    // }
                }
        // if(!loc_atp_index.count(ia1))
        //     continue;
        // if (box == Vector3_Order<int>({0, 0, 1}))continue;
        bool keep = loc_atp_index.count(ia1) && (*cs_ptr).absmax() >= threshold;
        if (keep)
            set_lri_coeff(ctx, keyword, ctx.opts.parallel_routing, ia1, ia2,
                                 n_i, n_j, n_mu, R, cs_ptr->c);
        // cout<<cs_ptr->nr<<cs_ptr->nc<<endl;
        if (!keep)
        {
            cs_discard++;
        }
    }
    infile.close();
    return cs_discard;
}

static size_t handle_Cs_file_binary(ReaderContext &ctx, const string &file_path, double threshold,
                                    const std::vector<atpair_t> &local_atpair,
                                    const string &keyword)
{
    using namespace std;

    librpa_int::require_readable_file(file_path);
    if (has_Cs_binary_v1_layout(file_path))
    {
        set<size_t> loc_atp_index;
        for (auto &lap : local_atpair)
        {
            loc_atp_index.insert(lap.first);
            loc_atp_index.insert(lap.second);
        }

        int natom = 0;
        int ncell = 0;
        std::streamoff file_size = 0;
        const auto records = read_Cs_binary_v1_header_or_throw(
            file_path, natom, ncell, file_size);

        const auto &basis_wfc = ctx.state.nbs_wfc;
        const auto &basis_aux = target_aux_basis_sizes_for_keyword(ctx, keyword);
        validate_Cs_binary_v1_blocks(
            file_path, natom, file_size, records, basis_wfc, basis_aux);

        ifstream infile(file_path, std::ios::in | std::ios::binary);
        if (!infile.good())
        {
            throw std::logic_error("Failed to open " + file_path);
        }

        size_t cs_discard = 0;
        for (const auto &block: records)
        {
            const int ia1 = block.ia1 - 1;
            const int ia2 = block.ia2 - 1;
            const bool keep =
                loc_atp_index.count(static_cast<size_t>(ia1)) && block.max_abs >= threshold;
            if (!keep)
            {
                cs_discard++;
                continue;
            }

            int n_i = 0;
            int n_j = 0;
            int n_mu = 0;
            get_Cs_binary_v1_block_dimensions(
                block, basis_wfc, basis_aux, file_path, n_i, n_j, n_mu);

            std::streamoff data_size = 0;
            if (!get_Cs_binary_data_size(n_i, n_j, n_mu, data_size))
            {
                throw std::runtime_error("Invalid Cs reader v1 block dimensions in: " +
                                         file_path);
            }

            shared_ptr<matrix> cs_ptr = make_shared<matrix>();
            cs_ptr->create(n_i * n_j, n_mu);
            infile.seekg(static_cast<std::streamoff>(block.offset), std::ios::beg);
            infile.read((char *) cs_ptr->c, data_size);
            if (!infile.good())
            {
                throw std::runtime_error("Failed to read Cs reader v1 payload in: " +
                                         file_path);
            }

            int R[3] = {block.R[0], block.R[1], block.R[2]};
            set_lri_coeff(ctx, keyword, ctx.opts.parallel_routing, ia1, ia2,
                                 n_i, n_j, n_mu, R, cs_ptr->c);
        }
        return cs_discard;
    }

    set<size_t> loc_atp_index;
    for (auto &lap : local_atpair)
    {
        loc_atp_index.insert(lap.first);
        loc_atp_index.insert(lap.second);
    }
    // cout<<"READING Cs from file: "<<file_path<<"  Cs_first_size:
    // "<<loc_atp_index.size()<<endl;
    // map<size_t,map<size_t,map<Vector3_Order<int>,std::shared_ptr<matrix>>>> Cs_m;
    size_t cs_discard = 0;
    ifstream infile;
    int dims[8];
    int n_apcell_file;
    int natom, ncell;

    infile.open(file_path, std::ios::in | std::ios::binary);
    infile.read((char *)&natom, sizeof(int));
    infile.read((char *)&ncell, sizeof(int));
    infile.read((char *)&n_apcell_file, sizeof(int));

    int R[3];

    for (int i = 0; i < n_apcell_file; i++)
    {
        infile.read((char *)&dims[0], 8 * sizeof(int));
        // cout<<ic_1<<mu_s<<endl;
        int ia1 = dims[0] - 1;
        int ia2 = dims[1] - 1;
        R[0] = dims[2];
        R[1] = dims[3];
        R[2] = dims[4];
        int n_i = dims[5];
        int n_j = dims[6];
        int n_mu = dims[7];

        // cout<< ia1<<ia2<<box<<endl;

        shared_ptr<matrix> cs_ptr = make_shared<matrix>();
        cs_ptr->create(n_i * n_j, n_mu);
        infile.read((char *)cs_ptr->c, n_i * n_j * n_mu * sizeof(double));
        bool keep = loc_atp_index.count(ia1) && (*cs_ptr).absmax() >= threshold;
        // cout << (*cs_ptr).absmax() << "\n";
        if (keep)
        {
            set_lri_coeff(ctx, keyword, ctx.opts.parallel_routing, ia1, ia2,
                                 n_i, n_j, n_mu, R, cs_ptr->c);
        }
        else
        {
            cs_discard++;
        }
    }
    return cs_discard;
}

size_t read_Cs(ReaderContext &ctx, const string &dir_path, double threshold,
               const std::vector<atpair_t> &local_atpair, const string keyword,
               int reader_version)
{
    using namespace std;

    if (reader_version < 0)
    {
        reader_version = detect_Cs_reader_version(ctx, dir_path, keyword);
    }
    librpa_int::global::lib_printf_root("LRI coefficient reader (%s): %s\n", keyword.c_str(),
                                        reader_version == 1 ? "binary v1" : "legacy");

    size_t cs_discard = 0;
    if (reader_version == 1)
    {
        const auto files = discover_Cs_files_for_keyword(ctx, dir_path, keyword);
        if (files.empty())
        {
            throw std::logic_error(
                "No LRI coefficient reader v1 files found with prefix " +
                keyword);
        }
        for (const auto &fn: files)
        {
            if (!has_Cs_binary_v1_layout(fn))
            {
                throw std::logic_error(
                    "LRI coefficient reader v1 expected a valid v1 header in: " + fn);
            }
            cs_discard += handle_Cs_file_binary(ctx, fn, threshold, local_atpair, keyword);
        }
        return cs_discard;
    }
    if (reader_version != 0)
    {
        throw_unknown_Cs_reader_version(reader_version);
    }

    // cout << "Begin to read Cs" << endl;
    // cout << "cs_threshold:  " << threshold << endl;
    bool binary;
    bool binary_checked = false;

    const auto files = discover_Cs_files_for_keyword(ctx, dir_path, keyword);
    if (files.empty())
    {
        throw std::logic_error("No LRI coefficient files found with prefix " + keyword);
    }
    for (const auto &fn: files)
    {
        if (has_Cs_binary_v1_layout(fn))
        {
            throw_Cs_v1_requires_reader_version(fn);
        }
        if (!binary_checked)
        {
            binary = check_Cs_file_binary(fn);
            binary_checked = true;
            if (ctx.comm.myid == 0)
            {
                if (binary)
                {
                    cout << "Unformatted binary Cs files detected" << endl;
                }
                else
                {
                    cout << "ASCII format Cs files detected" << endl;
                }
            }
        }
        if (binary)
        {
            cs_discard += handle_Cs_file_binary(ctx, fn, threshold, local_atpair, keyword);
        }
        else
        {
            cs_discard += handle_Cs_file(ctx, fn, threshold, local_atpair, keyword);
        }
    }
    // initialize basis set object
    // librpa_int::atomic_basis_wfc.set(atom_nw);
    // librpa_int::atomic_basis_abf.set(atom_mu);

    // atom_mu_part_range.resize(atom_mu.size());
    // atom_mu_part_range[0]=0;
    // for(int I=1;I!=atom_mu.size();I++)
    //     atom_mu_part_range[I]=atom_mu.at(I-1)+atom_mu_part_range[I-1];

    // N_all_mu=atom_mu_part_range[natom-1]+atom_mu[natom-1];
    // init_N_all_mu(); // FIXME: backward compat

    // for(int i=0;i!=atom_mu_part_range.size();i++)
    //     cout<<" atom_mu_part_range ,i: "<<i<<"    "<<atom_mu_part_range[i]<<endl;

    // cout << "Finish read Cs" << endl;
    return cs_discard;
}

std::vector<size_t> handle_Cs_file_dry(const string &file_path, double threshold)
{
    using namespace std;
    using namespace librpa_int;

    require_readable_file(file_path);
    std::vector<size_t> Cs_ids_keep;
    string natom_s, ncell_s, ia1_s, ia2_s, ic_1, ic_2, ic_3, i_s, j_s, mu_s, Cs_ele;
    ifstream infile;
    infile.open(file_path);
    infile >> natom_s >> ncell_s;
    // int natom = stoi(natom_s);
    // int ncell = stoi(ncell_s);

    size_t id = 0;
    // int R[3];

    while (infile.peek() != EOF)
    {
        infile >> ia1_s;
        if (infile.peek() == EOF) break;
        infile >> ia2_s >> ic_1 >> ic_2 >> ic_3 >> i_s >> j_s >> mu_s;
        int n_i = stoi(i_s);
        int n_j = stoi(j_s);
        int n_mu = stoi(mu_s);
        // int ia1 = stoi(ia1_s) - 1;
        // int ia2 = stoi(ia2_s) - 1;
        // R[0] = stoi(ic_1);
        // R[1] = stoi(ic_2);
        // R[2] = stoi(ic_3);
        // assign basis
        // set_ao_basis_aux(ia1, ia2, n_i, n_j, n_mu, R, nullptr, 1);

        double maxval = -1.0;
        for (int i = 0; i != n_i; i++)
            for (int j = 0; j != n_j; j++)
                for (int mu = 0; mu != n_mu; mu++)
                {
                    infile >> Cs_ele;
                    maxval = std::max(maxval, std::abs(stod(Cs_ele)));
                }
        // librpa_int::global::ofs_myid << id << " " << ia1 << " " << ia2 << " (" << ic_1 << ","
        //                              << ic_2 << "," << ic_3 << ") " << maxval << " keep? "
        //                              << (maxval >= threshold) << endl;
        if (maxval >= threshold) Cs_ids_keep.push_back(id);
        id++;
    }
    librpa_int::global::ofs_myid << file_path << ": " << Cs_ids_keep << endl;
    infile.close();
    return Cs_ids_keep;
}

std::vector<size_t> handle_Cs_file_binary_dry(const string &file_path, double threshold)
{
    librpa_int::require_readable_file(file_path);
    if (has_Cs_binary_v1_layout(file_path))
    {
        int natom = 0;
        int ncell = 0;
        std::streamoff file_size = 0;
        const auto records = read_Cs_binary_v1_header_or_throw(
            file_path, natom, ncell, file_size);

        std::vector<size_t> Cs_ids_keep;
        for (std::size_t i = 0; i != records.size(); ++i)
        {
            if (records[i].max_abs >= threshold) Cs_ids_keep.push_back(i);
        }
        return Cs_ids_keep;
    }

    std::vector<size_t> Cs_ids_keep;
    ifstream infile;
    int dims[8];
    int n_apcell_file;
    // int n_processed = 0;
    // int R[3];
    int natom, ncell;

    infile.open(file_path, std::ios::in | std::ios::binary);
    infile.read((char *)&natom, sizeof(int));
    infile.read((char *)&ncell, sizeof(int));
    infile.read((char *)&n_apcell_file, sizeof(int));

    for (int i_file = 0; i_file < n_apcell_file; i_file++)
    {
        infile.read((char *)&dims[0], 8 * sizeof(int));
        // cout<<ic_1<<mu_s<<endl;
        // const int ia1 = dims[0] - 1;
        // const int ia2 = dims[1] - 1;
        // R[0] = dims[2];
        // R[1] = dims[3];
        // R[2] = dims[4];
        const int n_i = dims[5];
        const int n_j = dims[6];
        const int n_mu = dims[7];
        // set_ao_basis_aux(ia1, ia2, n_i, n_j, n_mu, R, nullptr, 1);

        matrix mat(n_i * n_j, n_mu);
        infile.read((char *)mat.c, n_i * n_j * n_mu * sizeof(double));
        double maxval = mat.absmax();
        // n_processed++;
        if (maxval >= threshold)
        {
            Cs_ids_keep.push_back(i_file);
#ifdef LIBRPA_DEBUG
            // LIBRPA::envs::ofs_myid << i_file << " (" << ic1 << "," << ic2 << "," << ic3 << ")
            // "
            // << maxval << " kept, maxval: " << maxval << endl;
#endif
        }
    }
    // LIBRPA::envs::ofs_myid << file_path << ": kept " << Cs_ids_keep.size() << " of " <<
    // n_processed << endl;
#ifdef LIBRPA_DEBUG
    // librpa_int::envs::ofs_myid << Cs_ids_keep << endl;
#endif
    infile.close();
    return Cs_ids_keep;
}

static size_t handle_Cs_file_by_ids(ReaderContext &ctx, const std::string &file_path, double threshold,
                                    const std::vector<size_t> &ids, const std::string keyword)
{
    using namespace std;
    size_t cs_discard = 0;
    string natom_s, ncell_s, ia1_s, ia2_s, ic_1, ic_2, ic_3, i_s, j_s, mu_s, Cs_ele;
    librpa_int::require_readable_file(file_path);
    ifstream infile;
    // int natom, ncell;
    infile.open(file_path);
    infile >> natom_s >> ncell_s;
    // natom = stoi(natom_s);
    // ncell = stoi(ncell_s);
    /* cout<<"  Natom  Ncell  "<<natom<<"  "<<ncell<<endl; */
    // for(int loop=0;loop!=natom*natom*ncell;loop++)
    size_t id = 0;
    int R[3];

    while (infile.peek() != EOF)
    {
        infile >> ia1_s >> ia2_s >> ic_1 >> ic_2 >> ic_3 >> i_s;
        if (infile.peek() == EOF) break;
        // cout << " ia1_s,ia2_s: " << ia1_s << "  " << ia2_s << endl;
        infile >> j_s >> mu_s;
        // cout<<ic_1<<mu_s<<endl;
        int ia1 = stoi(ia1_s) - 1;
        int ia2 = stoi(ia2_s) - 1;
        R[0] = stoi(ic_1);
        R[1] = stoi(ic_2);
        R[2] = stoi(ic_3);
        int n_i = stoi(i_s);
        int n_j = stoi(j_s);
        int n_mu = stoi(mu_s);

        if (std::find(ids.cbegin(), ids.cend(), id) != ids.cend())
        {
            shared_ptr<matrix> cs_ptr = make_shared<matrix>();
            cs_ptr->create(n_i * n_j, n_mu);

            for (int i = 0; i != n_i; i++)
                for (int j = 0; j != n_j; j++)
                    for (int mu = 0; mu != n_mu; mu++)
                    {
                        infile >> Cs_ele;
                        (*cs_ptr)(i *n_j + j, mu) = stod(Cs_ele);
                    }
            set_lri_coeff(ctx, keyword, ctx.opts.parallel_routing, ia1, ia2,
                                 n_i, n_j, n_mu, R, cs_ptr->c);
        }
        else
        {
            // set_ao_basis_aux(ia1, ia2, n_i, n_j, n_mu, R, nullptr, 1, keyword);

            double maxval = -1.0;
            for (int i = 0; i != n_i; i++)
                for (int j = 0; j != n_j; j++)
                    for (int mu = 0; mu != n_mu; mu++)
                    {
                        infile >> Cs_ele;
                        maxval = std::max(maxval, std::abs(stod(Cs_ele)));
                    }
            if (maxval < threshold) cs_discard++;
        }
        id++;
    }
    infile.close();
    return cs_discard;
}


static size_t handle_Cs_file_binary_by_ids(ReaderContext &ctx, const string &file_path, double threshold,
                                           const std::vector<size_t> &ids,
                                           const string &keyword)
{
    using namespace std;

    librpa_int::require_readable_file(file_path);
    if (has_Cs_binary_v1_layout(file_path))
    {
        int natom = 0;
        int ncell = 0;
        std::streamoff file_size = 0;
        const auto records = read_Cs_binary_v1_header_or_throw(
            file_path, natom, ncell, file_size);

        const auto &basis_wfc = ctx.state.nbs_wfc;
        const auto &basis_aux = target_aux_basis_sizes_for_keyword(ctx, keyword);
        validate_Cs_binary_v1_blocks(
            file_path, natom, file_size, records, basis_wfc, basis_aux);

        ifstream infile(file_path, std::ios::in | std::ios::binary);
        if (!infile.good())
        {
            throw std::logic_error("Failed to open " + file_path);
        }

        size_t cs_discard = 0;
        for (std::size_t i = 0; i != records.size(); ++i)
        {
            const auto &block = records[i];
            const int ia1 = block.ia1 - 1;
            const int ia2 = block.ia2 - 1;

            if (std::find(ids.cbegin(), ids.cend(), i) == ids.cend())
            {
                cs_discard++;
                continue;
            }

            int n_i = 0;
            int n_j = 0;
            int n_mu = 0;
            get_Cs_binary_v1_block_dimensions(
                block, basis_wfc, basis_aux, file_path, n_i, n_j, n_mu);

            std::streamoff data_size = 0;
            if (!get_Cs_binary_data_size(n_i, n_j, n_mu, data_size))
            {
                throw std::runtime_error("Invalid Cs reader v1 block dimensions in: " +
                                         file_path);
            }

            shared_ptr<matrix> cs_ptr = make_shared<matrix>();
            cs_ptr->create(n_i * n_j, n_mu);
            infile.seekg(static_cast<std::streamoff>(block.offset), std::ios::beg);
            infile.read((char *) cs_ptr->c, data_size);
            if (!infile.good())
            {
                throw std::runtime_error("Failed to read Cs reader v1 payload in: " +
                                         file_path);
            }

            int R[3] = {block.R[0], block.R[1], block.R[2]};
            set_lri_coeff(ctx, keyword, ctx.opts.parallel_routing, ia1, ia2,
                                 n_i, n_j, n_mu, R, cs_ptr->c);
        }
        return cs_discard;
    }

    ifstream infile;
    int dims[8];
    int n_apcell_file;
    int natom, ncell;

    infile.open(file_path, std::ios::in | std::ios::binary);
    infile.read((char *)&natom, sizeof(int));
    infile.read((char *)&ncell, sizeof(int));
    infile.read((char *)&n_apcell_file, sizeof(int));
    size_t cs_discard = 0;

    int R[3];

    for (int i = 0; i < n_apcell_file; i++)
    {
        infile.read((char *)&dims[0], 8 * sizeof(int));
        // cout<<ic_1<<mu_s<<endl;
        int ia1 = dims[0] - 1;
        int ia2 = dims[1] - 1;
        R[0] = dims[2];
        R[1] = dims[3];
        R[2] = dims[4];
        int n_i = dims[5];
        int n_j = dims[6];
        int n_mu = dims[7];

        if (std::find(ids.cbegin(), ids.cend(), static_cast<size_t>(i)) != ids.cend())
        {
            shared_ptr<matrix> cs_ptr = make_shared<matrix>();
            cs_ptr->create(n_i * n_j, n_mu);
            infile.read((char *) cs_ptr->c, n_i * n_j * n_mu * sizeof(double));
            set_lri_coeff(ctx, keyword, ctx.opts.parallel_routing, ia1, ia2,
                                 n_i, n_j, n_mu, R, cs_ptr->c);
            // debug output
            // ofs_myid << "routing " << ctx.opts.parallel_routing << " ia1 " << ia1 << " ia2 "
            //          << ia2 << " R " << R[0] << " " << R[1] << " " << R[2] << " n_i " << n_i
            //          << " n_j " << n_j << " n_mu " << n_mu << endl;
            // print_matrix("cs_ptr->c", *cs_ptr, ofs_myid, true);
        }
        else
        {
            infile.seekg(n_i * n_j * n_mu * sizeof(double), ios::cur);
            cs_discard++;
        }
    }
    infile.close();
    return cs_discard;
}

size_t read_Cs_evenly_distribute(ReaderContext &ctx, const string &dir_path, double threshold, int myid, int nprocs,
                                 const string keyword, int reader_version)
{
    using namespace std;
    using namespace librpa_int;
    using namespace librpa_int::global;

    if (reader_version < 0)
    {
        reader_version = detect_Cs_reader_version(ctx, dir_path, keyword);
    }
    if (myid == 0)
    {
        lib_printf("LRI coefficient reader (%s): %s\n", keyword.c_str(),
                   reader_version == 1 ? "binary v1" : "legacy");
    }

    if (reader_version == 1)
    {
        auto files = discover_Cs_files_for_keyword(ctx, dir_path, keyword);
        if (files.empty())
        {
            throw std::logic_error(
                "No LRI coefficient reader v1 files found with prefix " + keyword);
        }

        size_t cs_discard = 0;
        profiler.start("handle_Cs_file_dry");
        const auto &basis_wfc = ctx.state.nbs_wfc;
        const auto &basis_aux = target_aux_basis_sizes_for_keyword(ctx, keyword);
        const auto tasks = make_Cs_binary_v1_read_tasks(files, threshold, basis_wfc, basis_aux);
        const auto tasks_this_proc = select_Cs_binary_v1_tasks_for_rank(tasks, myid, nprocs);
        cs_discard = tasks.size() - tasks_this_proc.size();
        ofs_myid << "Cs reader v1 scheduled " << tasks_this_proc.size()
                 << " of " << tasks.size() << " kept block(s) across "
                 << files.size() << " file(s)" << endl;
        profiler.stop("handle_Cs_file_dry");
        if (myid == 0) lib_printf("Finished Cs filtering\n");

        profiler.start("handle_Cs_file");
        read_Cs_binary_v1_tasks(ctx, tasks_this_proc, keyword);
        profiler.stop("handle_Cs_file");
        if (myid == 0) lib_printf("Finished Cs parsing\n");
        return cs_discard;
    }
    if (reader_version != 0)
    {
        throw_unknown_Cs_reader_version(reader_version);
    }

    size_t cs_discard = 0;
    auto files = discover_Cs_files_for_keyword(ctx, dir_path, keyword);
    unordered_map<string, std::vector<size_t>> files_Cs_ids;
    unordered_map<string, std::vector<size_t>> files_Cs_ids_this_proc;
    bool binary;
    bool binary_checked = false;

    profiler.start("handle_Cs_file_dry");
    if (files.empty())
    {
        profiler.stop("handle_Cs_file_dry");
        throw std::logic_error("No LRI coefficient files found with prefix " + keyword);
    }
    for (const auto &file_path: files)
    {
        if (has_Cs_binary_v1_layout(file_path))
        {
            profiler.stop("handle_Cs_file_dry");
            throw_Cs_v1_requires_reader_version(file_path);
        }
        if (!binary_checked)
        {
            binary = check_Cs_file_binary(file_path);
            binary_checked = true;
        }
    }

    const int nfiles = files.size();
    // cout << nfiles << "\n";

    // TODO: the IO can be improved, in two possible ways
    // 1. Each MPI task reads only a subset of files, instead of all files.
    // 2. Parallel reading for each file. This may be more efficient, but would be more difficult to
    // implement
    for (int i_fn = 0; i_fn != nfiles; i_fn++)
    {
        // Let each MPI process read different files at one time
        auto i_fn_myid = (i_fn + myid * nfiles / nprocs) % files.size();
        const auto &fn = files[i_fn_myid];
        ofs_myid << "Reading " << fn << " binary ? " << boolalpha << binary << endl;
        std::vector<size_t> ids_keep_this_file;
        if (binary)
        {
            ids_keep_this_file = handle_Cs_file_binary_dry(fn, threshold);
        }
        else
        {
            ids_keep_this_file = handle_Cs_file_dry(fn, threshold);
        }
        files_Cs_ids[fn] = ids_keep_this_file;
    }

    // Filter out the Cs to be actually read in each process
    int id_total = 0;
    for (int i_fn = 0; i_fn < nfiles; i_fn++)
    {
        const auto &fn = files[i_fn];
        const auto &ids_this_file = files_Cs_ids[fn];
        const int n_ids = ids_this_file.size();
        for (int id = 0; id != n_ids; id++)
        {
            if (id_total % nprocs == myid) files_Cs_ids_this_proc[fn].push_back(ids_this_file[id]);
            id_total++;
        }
    }
    profiler.stop("handle_Cs_file_dry");
    if (myid == 0) lib_printf("Finished Cs filtering\n");

    profiler.start("handle_Cs_file");
    // cout << files_Cs_ids_this_proc.size() << "\n";
    ofs_myid << "Number of Cs files to process: " << files_Cs_ids_this_proc.size() << "\n";
    for (const auto &fn_ids: files_Cs_ids_this_proc)
    {
        ofs_myid << fn_ids.first << " " << fn_ids.second << endl;
        if (binary)
        {
            cs_discard += handle_Cs_file_binary_by_ids(ctx, fn_ids.first, threshold, fn_ids.second,
                                                       keyword);
        }
        else
        {
            cs_discard += handle_Cs_file_by_ids(ctx, fn_ids.first, threshold, fn_ids.second, keyword);
        }
    }
    profiler.stop("handle_Cs_file");

    // initialize basis set object
    // librpa_int::atomic_basis_wfc.set(atom_nw);
    // librpa_int::atomic_basis_abf.set(atom_mu);

    // atom_mu_part_range.resize(atom_mu.size());
    // atom_mu_part_range[0]=0;
    // for(int I=1;I!=atom_mu.size();I++)
    //     atom_mu_part_range[I]=atom_mu.at(I-1)+atom_mu_part_range[I-1];
    //
    // N_all_mu=atom_mu_part_range[natom-1]+atom_mu[natom-1];
    // cout << "Done\n";
    if (myid == 0) lib_printf("Finished Cs parsing\n");

    return cs_discard;
}

void get_natom_ncell_from_first_Cs_file(ReaderContext &ctx, int &n_atom, int &n_cell, const string &dir_path)
{
    using namespace std;

    // cout<<file_path<<endl;
    ifstream infile;
    bool binary;

    string file_path = "";

    const auto files = discover_Cs_files_for_keyword(ctx, dir_path,
                                                     ctx.params.prefix_lri_coeff);
    if (!files.empty())
    {
        file_path = files.front();
    }
    if (file_path == "")
        throw std::runtime_error("Cs_data file is not found under dir_path: " + dir_path);

    binary = check_Cs_file_binary(file_path);
    if (ctx.comm.myid == 0)
    {
        if (binary)
        {
            cout << "Unformatted binary Cs files detected" << endl;
        }
        else
        {
            cout << "ASCII format Cs files detected" << endl;
        }
    }

    if (binary)
    {
        infile.open(file_path, std::ios::in | std::ios::binary);
        infile.read((char *)&n_atom, sizeof(int));
        infile.read((char *)&n_cell, sizeof(int));
        infile.close();
    }
    else
    {
        string natom_s, ncell_s;
        infile.open(file_path);
        infile >> natom_s >> ncell_s;
        // cout<<"  natom_s:"<<natom_s<<"  ncell_s: "<<ncell_s<<endl;
        n_atom = stoi(natom_s);
        n_cell = stoi(ncell_s);
        infile.close();
    }
}

static void get_basis_from_Cs(const string &file_path, std::map<int, size_t> &map_at_wfc, std::map<int, size_t> &map_at_aux)
{
    string natom_s, ncell_s, ia1_s, ia2_s, ic_1, ic_2, ic_3, i_s, j_s, mu_s, Cs_ele;
    librpa_int::require_readable_file(file_path);
    ifstream infile;
    infile.open(file_path);
    infile >> natom_s >> ncell_s;

    while (infile.peek() != EOF)
    {
        infile >> ia1_s >> ia2_s >> ic_1 >> ic_2 >> ic_3 >> i_s;
        if (infile.peek() == EOF)
            break;
        // cout << " ia1_s,ia2_s: " << ia1_s << "  " << ia2_s << endl;
        infile >> j_s >> mu_s;
        // cout<<ic_1<<mu_s<<endl;
        int ia1 = stoi(ia1_s) - 1;
        int ia2 = stoi(ia2_s) - 1;
        size_t n_i = stoi(i_s);
        size_t n_j = stoi(j_s);
        size_t n_mu = stoi(mu_s);
        map_at_wfc[ia1] = n_i;
        map_at_wfc[ia2] = n_j;
        map_at_aux[ia1] = n_mu;

        for (size_t i = 0; i != n_i; i++)
            for (size_t j = 0; j != n_j; j++)
                for (size_t mu = 0; mu != n_mu; mu++)
                    infile >> Cs_ele;
    }
    infile.close();
}

static void get_basis_from_Cs_binary(const string &file_path, std::map<int, size_t> &map_at_wfc, std::map<int, size_t> &map_at_aux)
{
    librpa_int::require_readable_file(file_path);
    if (has_Cs_binary_v1_layout(file_path))
    {
        throw std::runtime_error(
            "Cannot infer basis dimensions from Cs reader v1 file " + file_path +
            "; provide basis_out because v1 Cs block headers store offsets and max values, "
            "not n_i/n_j/n_mu");
    }

    ifstream infile;
    int dims[8];
    int n_apcell_file;
    int natom, ncell;
    const auto header_size = CS_BINARY_LEGACY_HEADER_SIZE;
    const auto block_header_size = CS_BINARY_LEGACY_BLOCK_HEADER_SIZE;

    infile.open(file_path, std::ios::in | std::ios::binary);
    if (!infile.good())
    {
        throw std::logic_error("Failed to open " + file_path);
    }
    infile.seekg(0, std::ios::end);
    const auto end_pos = infile.tellg();
    if (end_pos == std::streampos(-1))
        throw std::runtime_error("Failed to determine size of binary Cs file: " + file_path);
    const auto file_size = static_cast<std::streamoff>(end_pos);
    if (file_size < header_size)
        throw std::runtime_error("Binary Cs file is too small: " + file_path);

    infile.seekg(0, std::ios::beg);
    infile.read((char *) &natom, sizeof(int));
    infile.read((char *) &ncell, sizeof(int));
    infile.read((char *) &n_apcell_file, sizeof(int));
    if (!infile.good() || natom <= 0 || ncell < 0 || n_apcell_file < 0 ||
        n_apcell_file > (file_size - header_size) / block_header_size)
    {
        throw std::runtime_error("Invalid binary Cs header in: " + file_path);
    }

    std::streamoff pos = header_size;

    for (int i = 0; i < n_apcell_file; i++)
    {
        if (file_size - pos < block_header_size)
            throw std::runtime_error("Truncated binary Cs block header in: " + file_path);

        infile.seekg(pos, std::ios::beg);
        infile.read((char *) &dims[0], 8 * sizeof(int));
        if (!infile.good())
            throw std::runtime_error("Failed to read binary Cs block header in: " + file_path);
        pos += block_header_size;

        int ia1 = dims[0] - 1;
        int ia2 = dims[1] - 1;
        if (ia1 < 0 || ia1 >= natom || ia2 < 0 || ia2 >= natom)
            throw std::runtime_error("Invalid atom index in binary Cs block: " + file_path);

        std::streamoff data_size = 0;
        if (!get_Cs_binary_data_size(dims[5], dims[6], dims[7], data_size) ||
            file_size - pos < data_size)
        {
            throw std::runtime_error("Invalid binary Cs block dimensions in: " + file_path);
        }

        size_t n_i = static_cast<size_t>(dims[5]);
        size_t n_j = static_cast<size_t>(dims[6]);
        size_t n_mu = static_cast<size_t>(dims[7]);
        map_at_wfc[ia1] = n_i;
        map_at_wfc[ia2] = n_j;
        map_at_aux[ia1] = n_mu;
        pos += data_size;
    }
}


static void collect_basis_from_Cs_prefix(ReaderContext &ctx, const string &dir_path, const string &keyword,
                                         std::vector<size_t> &nbs_wfc,
                                         std::vector<size_t> &nbs_aux)
{
    using namespace librpa_int;

    bool binary;
    bool binary_checked = false;

    std::map<int, size_t> map_at_wfc;
    std::map<int, size_t> map_at_aux;

    const auto files = discover_Cs_files_for_keyword(ctx, dir_path, keyword);
    for (const auto &fn: files)
    {
        if (!binary_checked)
        {
            binary = check_Cs_file_binary(fn);
            binary_checked = true;
        }
        if (binary)
        {
            get_basis_from_Cs_binary(fn, map_at_wfc, map_at_aux);
        }
        else
        {
            get_basis_from_Cs(fn, map_at_wfc, map_at_aux);
        }
    }

    if (map_at_wfc.empty() || map_at_aux.empty())
        throw std::runtime_error("Cannot infer basis dimensions from Cs prefix " + keyword +
                                 " under: " + dir_path);
    if (map_at_wfc.size() != map_at_aux.size())
        throw std::runtime_error("Inconsistent wave-function and auxiliary basis dimensions from Cs prefix " +
                                 keyword);

    const int n_atoms = static_cast<int>(map_at_wfc.size());
    nbs_wfc.resize(n_atoms);
    nbs_aux.resize(n_atoms);
    for (int ia = 0; ia != n_atoms; ++ia)
    {
        nbs_wfc[ia] = map_at_wfc.at(ia);
        nbs_aux[ia] = map_at_aux.at(ia);
    }
}

std::vector<size_t> read_aux_basis_from_Cs(ReaderContext &ctx, const string &dir_path, const string &keyword)
{
    using namespace librpa_int;

    std::vector<size_t> nbs_wfc;
    std::vector<size_t> nbs_aux;
    int n_atoms = 0;

    if (ctx.comm.myid == 0)
    {
        collect_basis_from_Cs_prefix(ctx, dir_path, keyword, nbs_wfc, nbs_aux);
        n_atoms = static_cast<int>(nbs_aux.size());
    }

    ctx.comm.bcast(&n_atoms, 1, 0);
    if (ctx.comm.myid != 0)
    {
        nbs_aux.resize(n_atoms);
    }
    ctx.comm.bcast(nbs_aux.data(), n_atoms, 0);
    if (target_shrink_aux_basis_for_keyword(ctx, keyword))
        ctx.state.nbs_aux_shrink = nbs_aux;
    else
        ctx.state.nbs_aux = nbs_aux;
    return nbs_aux;
}

void read_basis_from_Cs(ReaderContext &ctx, const string &dir_path)
{
    using namespace librpa_int;

    global::lib_printf_root("Fallback reading basis information from Cs files under: %s\n", dir_path.c_str());

    int n_atoms;
    std::vector<size_t> nbs_wfc;
    std::vector<size_t> nbs_aux;

    // Let the master process reads and then broadcasts to others
    if (ctx.comm.myid == 0)
    {
        collect_basis_from_Cs_prefix(ctx, dir_path, ctx.params.prefix_lri_coeff,
                                     nbs_wfc, nbs_aux);
        n_atoms = static_cast<int>(nbs_wfc.size());
    }

    // Broadcast
    ctx.comm.bcast(&n_atoms, 1, 0);
    if (ctx.comm.myid != 0)
    {
        nbs_wfc.resize(n_atoms);
        nbs_aux.resize(n_atoms);
    }
    ctx.comm.bcast(nbs_wfc.data(), n_atoms, 0);
    ctx.comm.bcast(nbs_aux.data(), n_atoms, 0);
    ctx.h.set_ao_basis_wfc(nbs_wfc);
    ctx.h.set_ao_basis_aux(nbs_aux);
    ctx.state.nbs_wfc = nbs_wfc;
    ctx.state.nbs_aux = nbs_aux;
}

} // namespace librpa::reader
