#include "reader_basis.h"

#include "reader_context.h"

#include <librpa.hpp>
#include <librpa_enums.h>

#include "../../src/io/fs.h"
#include "../../src/io/global_io.h"
#include "../../src/utils/error.h"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace librpa::reader
{

using std::ifstream;
using std::string;

namespace
{

enum class BasisKind
{
    Wfc,
    Aux,
};

struct BasisHeader
{
    int ntypes = 0;
    bool split = false;
    std::string convention;
};

struct BasisLShells
{
    std::vector<std::vector<int>> wfc;
    std::vector<std::vector<int>> aux;
};

struct BasisInfo
{
    std::vector<size_t> nbs;
    std::vector<std::vector<int>> l_shells;
};

std::string normalize_basis_convention(const std::string &convention)
{
    std::string normalized = convention;
    std::transform(normalized.begin(), normalized.end(), normalized.begin(),
                   [](unsigned char ch) { return static_cast<char>(std::tolower(ch)); });
    normalized.erase(
        std::remove_if(normalized.begin(), normalized.end(),
                       [](unsigned char ch) {
                           return std::isspace(ch) != 0 || ch == '-' || ch == '_';
                       }),
        normalized.end());
    return normalized;
}

bool is_unset_basis_convention(const std::string &convention)
{
    return convention.empty() || convention == "unset" || convention == "unknown" ||
           convention == "fallback" || convention == "none";
}

std::string canonical_basis_convention_label(const std::string &convention)
{
    if (is_unset_basis_convention(convention))
        return "fallback";
    if (convention == "fhiaims")
        return "aims";
    return convention;
}

void parse_basis_convention(ReaderContext &ctx, const std::string &convention)
{
    const auto normalized = normalize_basis_convention(convention);
    const auto label = canonical_basis_convention_label(normalized);

    if (ctx.state.is_basis_convention_read)
    {
        if (label != ctx.state.basis_convention_label)
        {
            throw std::runtime_error("Inconsistent angular basis convention: " + convention +
                                     " (previously read " + ctx.state.basis_convention_label + ")");
        }
        return;
    }
    ctx.state.is_basis_convention_read = true;
    ctx.state.basis_convention_label = label;

    if (label == "fallback")
        return;

    bool known_convention = false;
    int bloch_phase, bloch_ratom;
    LibrpaAngularOrder order;
    LibrpaRshCoeff coeff_m_nega, coeff_m_posi;

    if (label == "aims")
    {
        known_convention = true;
        bloch_phase = -1;
        bloch_ratom = 0;
        order = LIBRPA_ANGULAR_ORDER_NATURAL;
        coeff_m_nega = LIBRPA_RSH_COEFF_1_M;
        coeff_m_posi = LIBRPA_RSH_COEFF_1_M;
    }
    if (label == "abacus")
    {
        known_convention = true;
        bloch_phase = -1;
        bloch_ratom = 0;
        order = LIBRPA_ANGULAR_ORDER_ABS_PM;
        coeff_m_nega = LIBRPA_RSH_COEFF_M_1;
        coeff_m_posi = LIBRPA_RSH_COEFF_1_M;
    }
    if (label == "openmx")
    {
        known_convention = true;
        bloch_phase = 1;
        bloch_ratom = 0;
        order = LIBRPA_ANGULAR_ORDER_OPENMX;
        coeff_m_nega = LIBRPA_RSH_COEFF_1_M;
        coeff_m_posi = LIBRPA_RSH_COEFF_M_1;
    }
    if (label == "pyscf")
    {
        known_convention = true;
        bloch_phase = 1;
        bloch_ratom = 0;
        order = LIBRPA_ANGULAR_ORDER_PYSCF;
        coeff_m_nega = LIBRPA_RSH_COEFF_1_M;
        coeff_m_posi = LIBRPA_RSH_COEFF_M_1;
    }
    if (known_convention)
    {
        ctx.h.set_basis_convention(bloch_phase, bloch_ratom, order,
                                       coeff_m_nega, coeff_m_posi);
        return;
    }
    throw std::runtime_error("Unknown angular basis convention: " + convention);
}

size_t parse_size_token(const std::string &token, const std::string &file_path)
{
    try
    {
        return static_cast<size_t>(std::stoull(token));
    }
    catch (const std::exception &)
    {
        throw LIBRPA_RUNTIME_ERROR("Invalid basis information header in " + file_path);
    }
}

int parse_int_token(const std::string &token, const std::string &file_path)
{
    try
    {
        return std::stoi(token);
    }
    catch (const std::exception &)
    {
        throw LIBRPA_RUNTIME_ERROR("Invalid basis information header in " + file_path);
    }
}

BasisHeader read_basis_header(std::istream &infile, const std::string &file_path)
{
    std::string line;
    while (std::getline(infile, line))
    {
        std::istringstream iss(line);
        std::vector<std::string> tokens;
        for (std::string token; iss >> token;)
        {
            tokens.push_back(token);
        }
        if (tokens.empty())
        {
            continue;
        }
        if (tokens.size() != 3 && tokens.size() != 4)
        {
            throw LIBRPA_RUNTIME_ERROR("Invalid basis information header in " + file_path);
        }

        BasisHeader header;
        header.ntypes = parse_int_token(tokens[0], file_path);
        parse_size_token(tokens[1], file_path);
        if (tokens.size() == 3)
        {
            header.split = true;
            header.convention = tokens[2];
        }
        else
        {
            parse_size_token(tokens[2], file_path);
            header.convention = tokens[3];
        }
        if (header.ntypes <= 0)
        {
            throw LIBRPA_RUNTIME_ERROR("Invalid basis information header in " + file_path);
        }
        return header;
    }
    throw LIBRPA_RUNTIME_ERROR("Invalid basis information header in " + file_path);
}

std::vector<int> read_basis_l_shells(
    std::istream &infile,
    const int nshell,
    const std::string &label,
    const std::string &file_path)
{
    if (nshell < 0)
    {
        throw LIBRPA_RUNTIME_ERROR("Invalid " + label + " shell-layout header in " + file_path);
    }

    std::vector<int> l_shells;
    l_shells.reserve(static_cast<std::size_t>(nshell));
    for (int ishell = 0; ishell < nshell; ++ishell)
    {
        int l = -1;
        if (!(infile >> l) || l < 0)
        {
            throw LIBRPA_RUNTIME_ERROR("Invalid " + label + " shell layout in " + file_path);
        }
        l_shells.push_back(l);
    }
    return l_shells;
}

std::vector<std::vector<int>> read_basis_type_l_shells(
    std::istream &infile,
    const int ntypes,
    const std::string &label,
    const std::string &file_path,
    int first_type,
    int first_nshell)
{
    std::vector<std::vector<int>> type_l_shells(static_cast<std::size_t>(ntypes));
    std::vector<bool> seen(static_cast<std::size_t>(ntypes), false);
    for (int itype = 0; itype < ntypes; ++itype)
    {
        int type = first_type;
        int nshell = first_nshell;
        if (itype != 0 && !(infile >> type >> nshell))
        {
            throw LIBRPA_RUNTIME_ERROR("Missing " + label + " shell-layout header in " + file_path);
        }
        const int type_index = type - 1;
        if (type_index < 0 || type_index >= ntypes || seen[static_cast<std::size_t>(type_index)])
        {
            throw LIBRPA_RUNTIME_ERROR("Invalid " + label + " shell-layout type in " + file_path);
        }
        type_l_shells[static_cast<std::size_t>(type_index)] =
            read_basis_l_shells(infile, nshell, label, file_path);
        seen[static_cast<std::size_t>(type_index)] = true;
    }
    return type_l_shells;
}

std::vector<std::vector<int>> assign_type_l_shells_to_atoms(ReaderContext &ctx,
    const std::vector<std::vector<int>> &type_l_shells,
    const int ntypes,
    const std::string &file_path)
{
    std::vector<std::vector<int>> shells;
    shells.reserve(ctx.state.atom_types.size());
    for (const int atom_type : ctx.state.atom_types)
    {
        if (atom_type < 0 || atom_type >= ntypes)
        {
            throw LIBRPA_RUNTIME_ERROR("Invalid atom type while assigning shell layouts in " + file_path);
        }
        shells.push_back(type_l_shells[static_cast<std::size_t>(atom_type)]);
    }
    return shells;
}

std::vector<std::vector<int>> read_single_basis_l_shells_tail(ReaderContext &ctx,
    std::istream &infile,
    const int ntypes,
    const std::string &label,
    const std::string &file_path)
{
    int type = 0;
    int nshell = 0;
    if (!(infile >> type >> nshell))
    {
        if (!infile.eof())
        {
            throw LIBRPA_RUNTIME_ERROR("Invalid shell-layout tail in " + file_path);
        }
        return {};
    }
    return assign_type_l_shells_to_atoms(ctx,
        read_basis_type_l_shells(infile, ntypes, label, file_path, type, nshell),
        ntypes,
        file_path);
}

BasisLShells read_combined_basis_l_shells_tail(ReaderContext &ctx,
    std::istream &infile,
    const int ntypes,
    const std::string &file_path)
{
    int type = 0;
    int nshell = 0;
    if (!(infile >> type >> nshell))
    {
        if (!infile.eof())
        {
            throw LIBRPA_RUNTIME_ERROR("Invalid shell-layout tail in " + file_path);
        }
        return {};
    }

    auto wfc_type_l_shells =
        read_basis_type_l_shells(infile, ntypes, "AO", file_path, type, nshell);
    if (!(infile >> type >> nshell))
    {
        throw LIBRPA_RUNTIME_ERROR("Missing ABF shell-layout header in " + file_path);
    }
    auto aux_type_l_shells =
        read_basis_type_l_shells(infile, ntypes, "ABF", file_path, type, nshell);

    return {assign_type_l_shells_to_atoms(ctx, wfc_type_l_shells, ntypes, file_path),
            assign_type_l_shells_to_atoms(ctx, aux_type_l_shells, ntypes, file_path)};
}

std::vector<size_t> read_basis_type_sizes(
    std::istream &infile,
    const BasisHeader &header,
    const BasisKind kind,
    const std::string &file_path)
{
    std::vector<size_t> type_nbs(static_cast<std::size_t>(header.ntypes), 0);
    std::vector<bool> seen(static_cast<std::size_t>(header.ntypes), false);
    for (int itype = 0; itype < header.ntypes; itype++)
    {
        int type = 0;
        size_t n_wfc = 0;
        size_t n_aux = 0;
        if (header.split)
        {
            if (!(infile >> type >> n_wfc))
            {
                throw LIBRPA_RUNTIME_ERROR("Invalid basis information body in " + file_path);
            }
            n_aux = n_wfc;
        }
        else if (!(infile >> type >> n_wfc >> n_aux))
        {
            throw LIBRPA_RUNTIME_ERROR("Invalid basis information body in " + file_path);
        }

        const int type_index = type - 1;
        if (type_index < 0 || type_index >= header.ntypes ||
            seen[static_cast<std::size_t>(type_index)])
        {
            throw LIBRPA_RUNTIME_ERROR("Invalid basis information type in " + file_path);
        }
        type_nbs[static_cast<std::size_t>(type_index)] =
            kind == BasisKind::Wfc ? n_wfc : n_aux;
        seen[static_cast<std::size_t>(type_index)] = true;
    }
    return type_nbs;
}

BasisInfo read_basis_info_from_basis_file(ReaderContext &ctx,
    const std::string &file_path,
    const BasisKind kind,
    const std::string &label)
{
    librpa_int::require_readable_file(file_path);
    ifstream infile(file_path);
    if (!infile.good())
    {
        throw LIBRPA_RUNTIME_ERROR("Failed to open basis information file " + file_path);
    }

    const int n_atoms = static_cast<int>(ctx.state.n_atoms);
    if (static_cast<size_t>(n_atoms) != ctx.state.atom_types.size())
    {
        throw LIBRPA_RUNTIME_ERROR("Number of atoms not consistent with the geometry file!");
    }

    const auto header = read_basis_header(infile, file_path);
    parse_basis_convention(ctx, header.convention);
    const auto type_nbs = read_basis_type_sizes(infile, header, kind, file_path);

    BasisInfo info;
    info.nbs.resize(static_cast<std::size_t>(n_atoms));
    for (int iat = 0; iat < n_atoms; iat++)
    {
        const auto type = ctx.state.atom_types[iat];
        if (type < 0 || type >= header.ntypes)
        {
            throw LIBRPA_RUNTIME_ERROR("Invalid atom type while assigning basis sizes in " + file_path);
        }
        info.nbs[static_cast<std::size_t>(iat)] = type_nbs[static_cast<std::size_t>(type)];
    }

    if (header.split)
    {
        info.l_shells = read_single_basis_l_shells_tail(ctx, infile, header.ntypes, label, file_path);
    }
    else
    {
        const auto shells = read_combined_basis_l_shells_tail(ctx, infile, header.ntypes, file_path);
        info.l_shells = kind == BasisKind::Wfc ? shells.wfc : shells.aux;
    }
    return info;
}

} // namespace

void reader_basis_wfc(ReaderContext &ctx, const std::string &file_path)
{
    librpa_int::global::lib_printf_root("Reading wave-function basis information file: %s\n",
                                        file_path.c_str());
    const auto basis_info = read_basis_info_from_basis_file(ctx, file_path, BasisKind::Wfc, "AO");
    ctx.h.set_ao_basis_wfc(basis_info.nbs, basis_info.l_shells);
    ctx.state.nbs_wfc = basis_info.nbs;
}

void reader_basis_aux(ReaderContext &ctx, const std::string &file_path)
{
    librpa_int::global::lib_printf_root("Reading auxiliary basis information file: %s\n",
                                        file_path.c_str());
    const auto basis_info = read_basis_info_from_basis_file(ctx, file_path, BasisKind::Aux, "ABF");
    ctx.h.set_ao_basis_aux(basis_info.nbs, basis_info.l_shells);
    ctx.state.nbs_aux = basis_info.nbs;
}

void reader_basis_aux_shrink(ReaderContext &ctx, const std::string &file_path)
{
    librpa_int::global::lib_printf_root("Reading shrink auxiliary basis information file: %s\n",
                                        file_path.c_str());
    const auto basis_info = read_basis_info_from_basis_file(ctx, file_path, BasisKind::Aux, "ABF");
    ctx.h.set_ao_basis_aux_shrink(basis_info.nbs, basis_info.l_shells);
    ctx.state.nbs_aux_shrink = basis_info.nbs;
}

void reader_basis(ReaderContext &ctx, const std::string &file_path)
{
    librpa_int::global::lib_printf_root("Reading basis information file: %s\n", file_path.c_str());
    const auto wfc_info = read_basis_info_from_basis_file(ctx, file_path, BasisKind::Wfc, "AO");
    const auto aux_info = read_basis_info_from_basis_file(ctx, file_path, BasisKind::Aux, "ABF");
    ctx.h.set_ao_basis_wfc(wfc_info.nbs, wfc_info.l_shells);
    ctx.h.set_ao_basis_aux(aux_info.nbs, aux_info.l_shells);
    ctx.state.nbs_wfc = wfc_info.nbs;
    ctx.state.nbs_aux = aux_info.nbs;
}

} // namespace librpa::reader
