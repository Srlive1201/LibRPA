#include "reader_structure.h"

#include "reader_context.h"

#include "../../src/io/fs.h"
#include "../../src/io/global_io.h"
#include "../../src/utils/error.h"

#include <algorithm>
#include <cctype>
#include <exception>
#include <fstream>
#include <vector>

namespace librpa::reader
{

namespace
{

std::string lowercase_token(std::string token)
{
    std::transform(token.begin(), token.end(), token.begin(),
                   [](unsigned char ch) { return static_cast<char>(std::tolower(ch)); });
    return token;
}

bool is_stru_symop_convention(const std::string &token)
{
    const auto convention = lowercase_token(token);
    return convention == "row" || convention == "col";
}

int parse_stru_int_token(const std::string &token, const std::string &context)
{
    try
    {
        std::size_t used = 0;
        const int value = std::stoi(token, &used);
        if (used == token.size())
        {
            return value;
        }
    }
    catch (const std::exception &)
    {
    }
    throw LIBRPA_RUNTIME_ERROR("Invalid integer in " + context + ": " + token);
}

double parse_stru_double_token(const std::string &token, const std::string &context)
{
    try
    {
        std::size_t used = 0;
        const double value = std::stod(token, &used);
        if (used == token.size())
        {
            return value;
        }
    }
    catch (const std::exception &)
    {
    }
    throw LIBRPA_RUNTIME_ERROR("Invalid floating-point value in " + context + ": " + token);
}

void require_stru_tail_tokens(const std::vector<std::string> &tokens,
                              const std::size_t pos,
                              const std::size_t count,
                              const std::string &context)
{
    if (pos > tokens.size() || tokens.size() - pos < count)
    {
        throw LIBRPA_RUNTIME_ERROR("Unexpected end of stru_out while reading " + context);
    }
}

bool is_stru_symop_header_at(const std::vector<std::string> &tokens, const std::size_t pos)
{
    return pos + 1 < tokens.size() && is_stru_symop_convention(tokens[pos + 1]);
}

bool is_stru_tail_boundary_at(const std::vector<std::string> &tokens, const std::size_t pos)
{
    return pos == tokens.size() || (pos < tokens.size() && is_stru_symop_header_at(tokens, pos));
}

std::size_t skip_legacy_stru_kpoint_section(ReaderContext &ctx, const std::vector<std::string> &tokens,
                                            std::size_t pos,
                                            const std::string &file_path)
{
    require_stru_tail_tokens(tokens, pos, 3, "legacy k-point grid");
    const int nk0 = parse_stru_int_token(tokens[pos], file_path);
    const int nk1 = parse_stru_int_token(tokens[pos + 1], file_path);
    const int nk2 = parse_stru_int_token(tokens[pos + 2], file_path);
    pos += 3;
    if (nk0 <= 0 || nk1 <= 0 || nk2 <= 0)
    {
        throw LIBRPA_RUNTIME_ERROR("Invalid legacy k-point grid in " + file_path);
    }

    const int nk_full = nk0 * nk1 * nk2;
    if (ctx.state.n_kpoints > 0 && ctx.state.n_kpoints <= nk_full)
    {
        const auto after_ibz_rows = pos + static_cast<std::size_t>(3 * ctx.state.n_kpoints);
        if (is_stru_tail_boundary_at(tokens, after_ibz_rows))
        {
            return after_ibz_rows;
        }
        const auto after_ibz_mapping = after_ibz_rows + static_cast<std::size_t>(nk_full);
        if (is_stru_tail_boundary_at(tokens, after_ibz_mapping))
        {
            return after_ibz_mapping;
        }
    }

    require_stru_tail_tokens(tokens, pos, static_cast<std::size_t>(3 * nk_full),
                             "legacy k-point rows");
    pos += static_cast<std::size_t>(3 * nk_full);
    if (is_stru_tail_boundary_at(tokens, pos))
    {
        return pos;
    }

    require_stru_tail_tokens(tokens, pos, static_cast<std::size_t>(nk_full),
                             "legacy k-point mapping");
    return pos + static_cast<std::size_t>(nk_full);
}

std::size_t read_stru_symops_from_tokens(const std::vector<std::string> &tokens,
                                         std::size_t pos,
                                         const std::string &file_path,
                                         int &n_symops,
                                         int &row_conv,
                                         std::vector<int> &rotmats,
                                         std::vector<double> &trans)
{
    require_stru_tail_tokens(tokens, pos, 2, "symmetry operation header");
    n_symops = parse_stru_int_token(tokens[pos], file_path);
    if (n_symops < 0)
    {
        throw LIBRPA_RUNTIME_ERROR("Invalid number of symmetry operations in " + file_path);
    }
    const auto convention = lowercase_token(tokens[pos + 1]);
    if (!is_stru_symop_convention(convention))
    {
        throw LIBRPA_RUNTIME_ERROR("Invalid symmetry operation convention in " + file_path
                                  + ": " + tokens[pos + 1]);
    }
    row_conv = convention == "row" ? 1 : 0;
    pos += 2;

    rotmats.clear();
    trans.clear();
    rotmats.reserve(static_cast<std::size_t>(9 * n_symops));
    trans.reserve(static_cast<std::size_t>(3 * n_symops));
    for (int isym = 0; isym != n_symops; ++isym)
    {
        require_stru_tail_tokens(tokens, pos, 12, "symmetry operation");
        for (int i = 0; i != 9; ++i)
        {
            rotmats.push_back(parse_stru_int_token(tokens[pos + static_cast<std::size_t>(i)],
                                                   file_path));
        }
        pos += 9;
        for (int i = 0; i != 3; ++i)
        {
            trans.push_back(parse_stru_double_token(tokens[pos + static_cast<std::size_t>(i)],
                                                    file_path));
        }
        pos += 3;
    }
    return pos;
}

void read_stru_tail_symops(ReaderContext &ctx, std::ifstream &infile, const std::string &file_path)
{
    std::vector<std::string> tokens;
    std::string token;
    while (infile >> token)
    {
        tokens.push_back(token);
    }
    if (tokens.empty())
    {
        return;
    }

    std::size_t pos = 0;
    if (tokens.size() < 2 || !is_stru_symop_convention(tokens[1]))
    {
        pos = skip_legacy_stru_kpoint_section(ctx, tokens, pos, file_path);
    }
    if (pos == tokens.size())
    {
        return;
    }
    if (pos + 1 >= tokens.size() || !is_stru_symop_convention(tokens[pos + 1]))
    {
        throw LIBRPA_RUNTIME_ERROR("Unexpected trailing data in " + file_path);
    }

    int n_symops = 0;
    int row_conv = 0;
    std::vector<int> rotmats;
    std::vector<double> trans;
    pos = read_stru_symops_from_tokens(tokens, pos, file_path, n_symops, row_conv,
                                       rotmats, trans);
    if (pos != tokens.size())
    {
        throw LIBRPA_RUNTIME_ERROR("Unexpected data after symmetry operations in " + file_path);
    }
    ctx.h.set_symmetry_operations(n_symops, row_conv,
                                      rotmats.empty() ? nullptr : rotmats.data(),
                                      trans.empty() ? nullptr : trans.data());
}

} // namespace

void reader_structure(ReaderContext &ctx, const std::string &file_path)
{
    using namespace librpa_int;
    global::lib_printf_root("Reading structure file: %s\n", file_path.c_str());

    require_readable_file(file_path);
    std::ifstream infile(file_path);
    if (!infile.good())
        throw LIBRPA_RUNTIME_ERROR("Fail to open structure file " + file_path);
    std::string x, y, z;

    std::vector<double> lat_mat(9);
    std::vector<double> G_mat(9);

    for (int i = 0; i < 3; i++)
    {
        infile >> x >> y >> z;
        lat_mat[i * 3] = stod(x);
        lat_mat[i * 3 + 1] = stod(y);
        lat_mat[i * 3 + 2] = stod(z);
    }

    for (int i = 0; i < 3; i++)
    {
        infile >> x >> y >> z;
        G_mat[i * 3] = stod(x);
        G_mat[i * 3 + 1] = stod(y);
        G_mat[i * 3 + 2] = stod(z);
    }

    ctx.h.set_latvec_and_G(lat_mat.data(), G_mat.data());

    infile >> ctx.state.n_atoms;
    const auto n_atoms = ctx.state.n_atoms;
    ctx.state.atom_types.resize(n_atoms);
    std::vector<double> coords(n_atoms * 3);
    int type;
    for (size_t iat = 0; iat < n_atoms; iat++)
    {
        for (int i = 0; i < 3; i++) infile >> coords[3 * iat + i];
        infile >> type;
        ctx.state.atom_types[iat] = type - 1;
    }
    ctx.h.set_atoms(ctx.state.atom_types, coords);
    read_stru_tail_symops(ctx, infile, file_path);
}

} // namespace librpa::reader
