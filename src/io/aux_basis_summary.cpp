#include "aux_basis_summary.h"

#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace librpa_int
{

namespace
{

std::string per_atom_text(const std::vector<std::size_t> &values)
{
    if (values.empty())
    {
        return "";
    }
    const auto first = values.front();
    for (const auto value : values)
    {
        if (value != first)
        {
            return "mixed";
        }
    }
    return std::to_string(first);
}

void append_separator(std::ostringstream &out)
{
    out << "+------+----------------+----------------+\n";
}

void append_row(std::ostringstream &out, const std::string &type,
                const std::string &large_per_atom, const std::string &small_per_atom)
{
    out << "| " << std::right << std::setw(4) << type << " | " << std::setw(14) << large_per_atom
        << " | " << std::setw(14) << small_per_atom << " |\n";
}

}  // namespace

std::string format_aux_basis_compression_summary(const std::vector<std::size_t> &large_nbs,
                                                 const std::vector<std::size_t> &small_nbs,
                                                 const std::map<atom_t, int> &atom_to_type)
{
    if (large_nbs.size() != small_nbs.size())
    {
        throw std::invalid_argument(
            "large and shrink auxiliary basis atom counts have different sizes");
    }

    std::map<int, std::vector<atom_t>> atoms_by_type;
    if (atom_to_type.size() == large_nbs.size())
    {
        for (const auto &[atom, type] : atom_to_type)
        {
            if (atom >= large_nbs.size())
            {
                throw std::invalid_argument(
                    "atom type map is inconsistent with auxiliary basis size");
            }
            atoms_by_type[type].push_back(atom);
        }
    }
    else
    {
        for (atom_t atom = 0; atom != large_nbs.size(); ++atom)
        {
            atoms_by_type[static_cast<int>(atom)].push_back(atom);
        }
    }

    std::ostringstream out;
    out << "Auxiliary basis compression summary:\n";
    append_separator(out);
    out << "| type | large ABF/atom | small ABF/atom |\n";
    append_separator(out);

    for (const auto &[type, atoms] : atoms_by_type)
    {
        std::vector<std::size_t> large_values;
        std::vector<std::size_t> small_values;
        large_values.reserve(atoms.size());
        small_values.reserve(atoms.size());
        for (const auto atom : atoms)
        {
            large_values.push_back(large_nbs[atom]);
            small_values.push_back(small_nbs[atom]);
        }
        append_row(out, std::to_string(type + 1), per_atom_text(large_values),
                   per_atom_text(small_values));
    }

    append_separator(out);
    return out.str();
}

}  // namespace librpa_int
