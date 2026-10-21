#include "crpa_window_input.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

#include "../math/rsh.h"

namespace librpa_int
{
namespace
{

using OrbitalComponent = std::pair<int, int>;  // (l, m), first radial shell

std::vector<OrbitalComponent> orbital_components(const std::string& orbitals)
{
    if (orbitals == "s") return {{0, 0}};
    if (orbitals == "p") return {{1, 0}, {1, 1}, {1, -1}};
    if (orbitals == "sp") return {{0, 0}, {1, 0}, {1, 1}, {1, -1}};
    if (orbitals == "d" || orbitals == "dp") return {{2, 0}, {2, 1}, {2, -1}, {2, 2}, {2, -2}};
    if (orbitals == "eg") return {{2, 0}, {2, 2}};
    if (orbitals == "t2g") return {{2, 1}, {2, -1}, {2, -2}};
    throw LIBRPA_INVALID_ARGUMENT("cRPA orbitals must be s, p, sp, d, dp, eg, or t2g");
}

std::string orbital_label(int l, int m)
{
    if (l == 0 && m == 0) return "s";
    if (l == 1)
    {
        if (m == 0) return "pz";
        return m == 1 ? "px" : "py";
    }
    switch (m)
    {
        case 0:
            return "dz2";
        case 1:
            return "dxz";
        case -1:
            return "dyz";
        case 2:
            return "dx2-y2";
        case -2:
            return "dxy";
    }
    throw LIBRPA_INVALID_ARGUMENT("invalid real-spherical orbital component");
}

int species_index(const std::vector<std::string>& labels, const std::string& selected)
{
    if (selected.empty() || std::count(labels.begin(), labels.end(), selected) != 1)
        throw LIBRPA_INVALID_ARGUMENT("cRPA species must identify one explicit species: " +
                                      selected);
    return static_cast<int>(std::find(labels.begin(), labels.end(), selected) - labels.begin());
}

int first_shell_offset(const std::vector<int>& shells, int l)
{
    int offset = 0;
    for (const int shell : shells)
    {
        if (shell == l) return offset;
        offset += 2 * shell + 1;
    }
    throw LIBRPA_INVALID_ARGUMENT("selected cRPA atom has no shell with l=" + std::to_string(l));
}

}  // namespace

std::vector<std::vector<std::vector<int>>> crpa_bands_in_energy_windows(
    const MeanField& meanfield, const std::vector<double>& edges_ha)
{
    if (!meanfield.initialized() || edges_ha.empty() || edges_ha.size() % 2 != 0)
        throw LIBRPA_INVALID_ARGUMENT(
            "cRPA requires a mean field and lower/upper energy-window pairs");
    for (std::size_t i = 0; i < edges_ha.size(); i += 2)
        if (!std::isfinite(edges_ha[i]) || !std::isfinite(edges_ha[i + 1]) ||
            edges_ha[i] > edges_ha[i + 1])
            throw LIBRPA_INVALID_ARGUMENT("invalid cRPA energy-window bounds");

    const auto& energies = meanfield.get_eigenvals();
    const int n_spins = meanfield.get_n_spins();
    const int nk = meanfield.get_n_kpoints();
    const int nb = meanfield.get_n_bands();
    if (energies.size() != static_cast<std::size_t>(n_spins))
        throw LIBRPA_INVALID_ARGUMENT("cRPA KS energy spin dimensions are incomplete");
    std::vector<std::vector<std::vector<int>>> bands(n_spins, std::vector<std::vector<int>>(nk));
    for (int spin = 0; spin < n_spins; ++spin)
    {
        if (energies[spin].nr != nk || energies[spin].nc != nb)
            throw LIBRPA_INVALID_ARGUMENT("cRPA KS energy dimensions are inconsistent");
        for (int k = 0; k < nk; ++k)
            for (int band = 0; band < nb; ++band)
            {
                const double energy = energies[spin](k, band);
                if (!std::isfinite(energy))
                    throw LIBRPA_INVALID_ARGUMENT("non-finite original KS energy");
                for (std::size_t edge = 0; edge < edges_ha.size(); edge += 2)
                    if (energy >= edges_ha[edge] && energy <= edges_ha[edge + 1])
                    {
                        bands[spin][k].push_back(band);
                        break;
                    }
            }
    }
    return bands;
}

namespace
{
std::vector<std::vector<std::vector<int>>> select_bands(const MeanField& meanfield,
                                                        const std::vector<double>& edges,
                                                        const std::vector<int>& explicit_bands)
{
    if (edges.empty() == explicit_bands.empty())
        throw LIBRPA_INVALID_ARGUMENT(
            "each cRPA subspace requires exactly one energy-window or explicit-band selection");
    if (!edges.empty()) return crpa_bands_in_energy_windows(meanfield, edges);
    auto bands = explicit_bands;
    std::sort(bands.begin(), bands.end());
    if (bands.front() < 0 || bands.back() >= meanfield.get_n_bands() ||
        std::adjacent_find(bands.begin(), bands.end()) != bands.end())
        throw LIBRPA_INVALID_ARGUMENT("cRPA explicit band indices must be distinct and in range");
    return std::vector<std::vector<std::vector<int>>>(
        meanfield.get_n_spins(), std::vector<std::vector<int>>(meanfield.get_n_kpoints(), bands));
}
}  // namespace

ComplexMatrix crpa_spinor_overlap(const ComplexMatrix& overlap, int n_spinor)
{
    if (overlap.nr <= 0 || overlap.nr != overlap.nc || (n_spinor != 1 && n_spinor != 2))
        throw LIBRPA_INVALID_ARGUMENT(
            "cRPA requires a square spatial overlap and one or two components");
    if (n_spinor == 1) return overlap;
    ComplexMatrix lifted(n_spinor * overlap.nr, n_spinor * overlap.nc);
    for (int s = 0; s < n_spinor; ++s)
        for (int i = 0; i < overlap.nr; ++i)
            for (int j = 0; j < overlap.nc; ++j)
                lifted(s * overlap.nr + i, s * overlap.nc + j) = overlap(i, j);
    return lifted;
}

CrpaWindowBuildResult build_crpa_window_input(
    const AtomicBasis& basis, const BasisConvention& convention, const std::vector<int>& atom_types,
    const std::vector<std::string>& species_labels, const std::vector<ComplexMatrix>& overlap_k,
    const std::vector<Vector3_Order<double>>& kfrac_list, const MeanField& meanfield,
    const CrpaWindowOptions& options)
{
    const int n_spinor = meanfield.get_n_spinor();
    const int nao = meanfield.get_n_aos();
    const int n_rows = nao * n_spinor;
    if (!basis.initialized() || !basis.has_l_shells() || !meanfield.initialized() ||
        (n_spinor != 1 && (n_spinor != 2 || meanfield.get_n_spins() != 1)) ||
        basis.nb_total != nao || atom_types.size() != basis.n_atoms || kfrac_list.empty() ||
        kfrac_list.size() != static_cast<std::size_t>(meanfield.get_n_kpoints()) ||
        overlap_k.size() != kfrac_list.size())
        throw LIBRPA_INVALID_ARGUMENT(
            "cRPA window construction requires matching scalar or two-component KS/basis/k data");
    // The downstream lattice Fourier transform uses the same cell-periodic convention.
    if (!is_basis_convention_set(convention) || convention.bloch_phase != -1 ||
        convention.bloch_ratom != 0)
        throw LIBRPA_INVALID_ARGUMENT(
            "cRPA window orbitals require a specified cell-periodic basis convention");

    const auto parent_m = orbital_components(options.parent_orbitals);
    const auto output_m = orbital_components(options.output_orbitals);
    const bool parent_p = options.parent_orbitals == "dp";
    const bool output_p = options.output_orbitals == "dp";
    if (output_p && !parent_p)
        throw LIBRPA_INVALID_ARGUMENT("dp output requires a dp parent orbital frame");
    for (const auto& m : output_m)
        if (std::find(parent_m.begin(), parent_m.end(), m) == parent_m.end())
            throw LIBRPA_INVALID_ARGUMENT("output orbitals are not a subset of the parent frame");

    const int correlated = species_index(species_labels, options.correlated_species);
    const int ligand =
        options.ligand_species.empty() ? -1 : species_index(species_labels, options.ligand_species);
    if ((parent_p && ligand < 0) || ligand == correlated)
        throw LIBRPA_INVALID_ARGUMENT("dp parent requires a distinct explicit ligand species");
    for (const int type : atom_types)
        if (type < 0 || type >= static_cast<int>(species_labels.size()))
            throw LIBRPA_INVALID_ARGUMENT("atom type is outside the species label table");

    CrpaWindowBuildResult result;
    if (options.use_orbital_weights)
    {
        if (!options.response_windows_ha.empty() || !options.response_bands.empty())
            throw LIBRPA_INVALID_ARGUMENT(
                "orbital-weight cRPA does not accept an A window or band list");
        result.response_bands.resize(meanfield.get_n_spins(),
                                     std::vector<std::vector<int>>(meanfield.get_n_kpoints()));
        result.response_weights.resize(meanfield.get_n_spins(),
                                       std::vector<std::vector<double>>(meanfield.get_n_kpoints()));
    }
    else
        result.response_bands =
            select_bands(meanfield, options.response_windows_ha, options.response_bands);
    result.orbital_bands =
        select_bands(meanfield, options.orbital_windows_ha, options.orbital_bands);
    std::vector<SiteOrbitalGroup> parent_sites;
    std::vector<int> trial_rows;
    std::vector<int> output_columns;
    int n_parent = 0;
    int n_output = 0;
    auto append_site = [&](int atom, const std::vector<OrbitalComponent>& components,
                           const std::vector<OrbitalComponent>& selected)
    {
        SiteOrbitalGroup parent;
        parent.label = species_labels[atom_types[atom]] + "-atom-" + std::to_string(atom);
        parent.atom_index = atom;
        parent.orb_start = n_parent;
        parent.n_orbitals = n_spinor * static_cast<int>(components.size());
        for (int s = 0; s < n_spinor; ++s)
            for (const auto& component : components)
            {
                const int l = component.first;
                const int m = component.second;
                const int shell = first_shell_offset(basis.get_l_shells(atom), l);
                const int local = shell + rsh_m_to_index(l, m, convention.order);
                if (local < 0 || local >= static_cast<int>(basis.get_atom_nb(atom)))
                    throw LIBRPA_INVALID_ARGUMENT(
                        "cRPA shell metadata exceeds the atom basis size");
                trial_rows.push_back(s * nao + static_cast<int>(basis.get_part_range()[atom]) +
                                     local);
                parent.orbital_labels.push_back(orbital_label(l, m) +
                                                (n_spinor == 2 ? (s == 0 ? ":up" : ":down") : ""));
            }
        parent_sites.push_back(parent);
        if (!selected.empty())
        {
            auto output = parent;
            output.orb_start = n_output;
            output.n_orbitals = n_spinor * static_cast<int>(selected.size());
            output.orbital_labels.clear();
            for (int s = 0; s < n_spinor; ++s)
                for (const auto& component : selected)
                {
                    const auto found = std::find(components.begin(), components.end(), component);
                    output_columns.push_back(n_parent + s * static_cast<int>(components.size()) +
                                             static_cast<int>(found - components.begin()));
                    output.orbital_labels.push_back(
                        orbital_label(component.first, component.second) +
                        (n_spinor == 2 ? (s == 0 ? ":up" : ":down") : ""));
                }
            result.sites.push_back(output);
            n_output += output.n_orbitals;
        }
        n_parent += parent.n_orbitals;
    };
    for (int atom = 0; atom < static_cast<int>(basis.n_atoms); ++atom)
        if (atom_types[atom] == correlated) append_site(atom, parent_m, output_m);
    if (parent_sites.empty())
        throw LIBRPA_INVALID_ARGUMENT("no atom matches the selected correlated species");
    if (parent_p)
    {
        const int before = n_parent;
        for (int atom = 0; atom < static_cast<int>(basis.n_atoms); ++atom)
            if (atom_types[atom] == ligand)
                append_site(atom, orbital_components("p"),
                            output_p ? orbital_components("p") : std::vector<OrbitalComponent>{});
        if (before == n_parent)
            throw LIBRPA_INVALID_ARGUMENT("no atom matches the selected ligand species");
    }

    ComplexMatrix trials(n_rows, n_parent);
    trials.zero_out();
    for (int column = 0; column < n_parent; ++column) trials(trial_rows[column], column) = 1.0;
    result.orbitals_spin_k.resize(meanfield.get_n_spins());
    result.parent_sites = parent_sites;
    result.output_columns = output_columns;
    result.parent_orbitals_spin_k.resize(meanfield.get_n_spins());
    const bool tracking = !options.previous_parent_spin_k.empty();
    if (!std::isfinite(options.tracking_min_capture) || options.tracking_min_capture <= 0 ||
        options.tracking_min_capture > 1 ||
        (tracking && options.previous_parent_spin_k.size() != meanfield.get_n_spins()))
        throw LIBRPA_INVALID_ARGUMENT("invalid cRPA parent-tracking dimensions or capture");
    result.parent_capture_min.assign(meanfield.get_n_spins(),
                                     std::vector<double>(meanfield.get_n_kpoints(), 1.0));
    for (int spin = 0; spin < meanfield.get_n_spins(); ++spin)
    {
        // No Fourier transform is needed until the API constructs the final output subspace.
        CorrelatedSubspace parent(
            parent_sites, kfrac_list, {{0, 0, 0}}, n_rows, meanfield.get_n_bands(), 1,
            options.gram_abs, options.gram_rel, options.gram_cond_max, options.s_abs, options.s_rel,
            options.s_cond_max, options.residual_tol, options.t_roundtrip_tol);
        for (int k = 0; k < meanfield.get_n_kpoints(); ++k)
        {
            ComplexMatrix wfc(meanfield.get_n_bands(), n_rows);
            for (int s = 0; s < n_spinor; ++s)
            {
                const auto* component = meanfield.find_wfc(spin, s, k);
                if (!component || component->nr != wfc.nr || component->nc != nao)
                    throw LIBRPA_INVALID_ARGUMENT(
                        "cRPA output construction requires all KS spinor components");
                for (int band = 0; band < wfc.nr; ++band)
                    for (int row = 0; row < nao; ++row)
                        wfc(band, s * nao + row) = (*component)(band, row);
            }
            if (overlap_k[k].nr != nao || overlap_k[k].nc != nao)
                throw LIBRPA_INVALID_ARGUMENT("cRPA spatial overlap has an inconsistent shape");
            const auto overlap = crpa_spinor_overlap(overlap_k[k], n_spinor);
            parent.set_S_k(k, overlap);
            const auto& reference =
                tracking ? options.previous_parent_spin_k.at(spin).at(k) : trials;
            if (reference.nr != n_rows || reference.nc != n_parent)
                throw LIBRPA_INVALID_ARGUMENT(
                    "cRPA previous parent must retain all parent columns");
            if (tracking)
            {
                const auto metric = transpose(reference, true) * overlap * reference;
                double residual = 0;
                for (int a = 0; a < n_parent; ++a)
                    for (int b = 0; b < n_parent; ++b)
                        residual += std::norm(metric(a, b) - (a == b ? 1.0 : 0.0));
                if (!std::isfinite(residual) ||
                    std::sqrt(residual) > options.residual_tol * std::sqrt(n_parent))
                    throw LIBRPA_INVALID_ARGUMENT("cRPA previous parent is not metric orthonormal");
            }
            parent.set_W_k(k, project_trial_orbitals_to_band_window(wfc, overlap, reference,
                                                                    result.orbital_bands[spin][k]));
            const auto phi = parent.build_phi(k);
            if (tracking)
            {
                const double capture = parent.get_kpt_diagnostics().at(k).gram_lambda_min;
                result.parent_capture_min[spin][k] = capture;
                if (capture < options.tracking_min_capture)
                    throw LIBRPA_INVALID_ARGUMENT(
                        "cRPA B loses a previous parent direction; widen the declared B "
                        "(minimum retained probability=" +
                        std::to_string(capture) + ")");
            }
            parent.build_spin_k(0, k, wfc);
            ComplexMatrix output(n_rows, n_output);
            for (int row = 0; row < n_rows; ++row)
                for (int column = 0; column < n_output; ++column)
                    output(row, column) = phi(row, output_columns[column]);
            if (options.use_orbital_weights)
            {
                const auto sc = overlap * transpose(wfc, false);
                const auto projector = transpose(output, true) * sc;
                auto& weights = result.response_weights[spin][k];
                weights.assign(meanfield.get_n_bands(), 0.0);
                const double tolerance = options.residual_tol * std::sqrt(n_rows);
                for (int band = 0; band < meanfield.get_n_bands(); ++band)
                {
                    std::complex<double> norm = 0;
                    for (int row = 0; row < n_rows; ++row)
                        norm += std::conj(wfc(band, row)) * sc(row, band);
                    if (!std::isfinite(norm.real()) || norm.real() <= 0 ||
                        std::abs(norm.imag()) > tolerance)
                        throw LIBRPA_INVALID_ARGUMENT("invalid KS metric norm in cRPA weights");
                    double weight = 0;
                    for (int a = 0; a < n_output; ++a) weight += std::norm(projector(a, band));
                    weight /= norm.real();
                    if (!std::isfinite(weight) || weight > 1.0 + tolerance)
                        throw LIBRPA_INVALID_ARGUMENT("cRPA projector probability exceeds one");
                    // Only remove accepted input/orthogonalization roundoff at one.
                    // Small positive weights are retained without a threshold.
                    weights[band] = std::min(1.0, weight);
                }
            }
            result.orbitals_spin_k[spin].push_back(std::move(output));
            result.parent_orbitals_spin_k[spin].push_back(phi);
        }
    }
    return result;
}

}  // namespace librpa_int
