#pragma once

// Historical upper-atom-block contraction controls, confined to this test.
// Production contracts full 2D BLACS kernels in crpa.cpp. The independent
// literal vertex/contraction oracles remain in the numerical test sources.
#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "../core/pair_density_vertex.h"

namespace crpa_test
{
using namespace librpa_int;
using PairInteractionAtomBlocks = std::map<int, std::map<int, ComplexMatrix>>;

inline bool finite(const std::complex<double> &value)
{
    return std::isfinite(value.real()) && std::isfinite(value.imag());
}

inline void validate_finite_matrix(const ComplexMatrix &m, const std::string &label)
{
    for (int i = 0; i != m.size; ++i)
    {
        if (!finite(m.c[i])) throw LIBRPA_RUNTIME_ERROR(label + " contains a non-finite value");
    }
}

inline void collective_stage_error(const MpiCommHandler &comm_h, const std::string &label,
                                   const std::string &local_error)
{
    int local_failed = local_error.empty() ? 0 : 1;
    int global_failed = 0;
    MPI_Allreduce(&local_failed, &global_failed, 1, MPI_INT, MPI_MAX, comm_h.comm);
    if (global_failed != 0)
        throw LIBRPA_RUNTIME_ERROR(
            label + " failed on at least one MPI rank: " +
            (local_error.empty() ? "another rank reported an error" : local_error));
}

inline void add_ordered_pair_block(ComplexMatrix &U, const ComplexMatrix &D_left,
                                   const ComplexMatrix &W, const ComplexMatrix &D_right,
                                   const int n_orbitals)
{
    const int npair = n_orbitals * n_orbitals;
    if (D_left.nr != npair || D_right.nr != npair || W.nr != D_left.nc || W.nc != D_right.nc)
        throw LIBRPA_INVALID_ARGUMENT(
            "ordered-pair contraction received inconsistent D/W block dimensions");

    U += contract_rectangular_ordered_pair_vertex(D_left, W, D_right, n_orbitals, n_orbitals);
}

inline ComplexMatrix contract_distributed_rectangular_ordered_pair_vertex(
    const PairVertexAtomBlocks &D_left, const PairVertexAtomBlocks &D_right,
    const PairInteractionAtomBlocks &W_local_upper, const int n_left, const int n_right,
    const MpiCommHandler &comm_h, const std::string &label)
{
    if (!comm_h.is_initialized())
        throw LIBRPA_INVALID_ARGUMENT(
            label + ": distributed rectangular contraction requires an initialized communicator");
    const int left_pairs = n_left * n_left;
    const int right_pairs = n_right * n_right;
    if (n_left <= 0 || n_right <= 0 || D_left.empty() || D_right.empty())
        throw LIBRPA_INVALID_ARGUMENT(label + ": invalid rectangular vertex dimensions or support");
    ComplexMatrix local(left_pairs, right_pairs);
    std::string local_error;
    for (const auto &[I, left_block] : D_left)
        for (const auto &[J, right_block] : D_right)
        {
            const bool forward = I <= J;
            const int row_atom = forward ? I : J;
            const int column_atom = forward ? J : I;
            const auto row_it = W_local_upper.find(row_atom);
            const bool local_owner = row_it != W_local_upper.end() &&
                                     row_it->second.find(column_atom) != row_it->second.end();
            int owners = 0;
            const int owner = local_owner ? 1 : 0;
            MPI_Allreduce(&owner, &owners, 1, MPI_INT, MPI_SUM, comm_h.comm);
            if (owners != 1)
            {
                local_error = label + ": each cross-site W block requires exactly one MPI owner";
                continue;
            }
            if (!local_owner) continue;
            const auto &stored = row_it->second.at(column_atom);
            const ComplexMatrix W = forward ? stored : transpose(stored, true);
            try
            {
                local += contract_rectangular_ordered_pair_vertex(left_block, W, right_block,
                                                                  n_left, n_right);
            }
            catch (const std::exception &error)
            {
                local_error = label + ": " + error.what();
            }
        }
    collective_stage_error(comm_h, label + " block contraction", local_error);
    if (local.size > std::numeric_limits<int>::max())
        throw LIBRPA_RUNTIME_ERROR(label + ": rectangular result exceeds MPI count capacity");
    MPI_Allreduce(MPI_IN_PLACE, local.c, local.size, MPI_CXX_DOUBLE_COMPLEX, MPI_SUM, comm_h.comm);
    validate_finite_matrix(local, label + " result");
    return local;
}

inline ComplexMatrix contract_ordered_pair_vertex(const PairVertexAtomBlocks &D_by_atom,
                                                  const PairInteractionAtomBlocks &W_upper,
                                                  const int n_orbitals)
{
    if (n_orbitals <= 0)
        throw LIBRPA_INVALID_ARGUMENT("ordered-pair contraction requires a positive orbital count");
    if (D_by_atom.empty() || W_upper.empty())
        throw LIBRPA_INVALID_ARGUMENT("ordered-pair contraction requires non-empty D and W blocks");
    const int npair = n_orbitals * n_orbitals;
    for (const auto &[I, D] : D_by_atom)
    {
        if (I < 0 || D.nr != npair || D.nc <= 0)
            throw LIBRPA_INVALID_ARGUMENT(
                "ordered-pair contraction received an invalid D atom block");
        validate_finite_matrix(D, "ordered-pair D block");
    }

    ComplexMatrix U(npair, npair);
    bool used_block = false;
    for (const auto &[I, J_W] : W_upper)
        for (const auto &[J, W] : J_W)
        {
            if (I > J)
                throw LIBRPA_INVALID_ARGUMENT(
                    "ordered-pair contraction expects upper-triangular W atom blocks");
            validate_finite_matrix(W, "ordered-pair W block");
            const auto it_I = D_by_atom.find(I);
            const auto it_J = D_by_atom.find(J);
            if (it_I == D_by_atom.end() || it_J == D_by_atom.end()) continue;
            if (W.nr != it_I->second.nc || W.nc != it_J->second.nc)
                throw LIBRPA_INVALID_ARGUMENT(
                    "ordered-pair contraction W dimensions do not match D blocks");
            if (I == J)
            {
                double scale = 1.0;
                double residual = 0.0;
                for (int i = 0; i != W.nr; ++i)
                    for (int j = 0; j != W.nc; ++j)
                    {
                        scale = std::max(scale, std::abs(W(i, j)));
                        residual = std::max(residual, std::abs(W(i, j) - std::conj(W(j, i))));
                    }
                if (residual > 1.0e-10 * scale)
                    throw LIBRPA_RUNTIME_ERROR(
                        "ordered-pair contraction diagonal W block is not Hermitian");
            }

            add_ordered_pair_block(U, it_I->second, W, it_J->second, n_orbitals);
            if (I != J)
            {
                const ComplexMatrix W_dagger = transpose(W, true);
                add_ordered_pair_block(U, it_J->second, W_dagger, it_I->second, n_orbitals);
            }
            used_block = true;
        }
    if (!used_block)
        throw LIBRPA_RUNTIME_ERROR("ordered-pair contraction found no W block on the D support");
    validate_finite_matrix(U, "ordered-pair U matrix");
    return U;
}

inline ComplexMatrix contract_distributed_ordered_pair_vertex(
    const PairVertexAtomBlocks &D_by_atom, const PairInteractionAtomBlocks &W_local_upper,
    const int n_orbitals, const MpiCommHandler &comm_h, const std::string &label)
{
    if (!comm_h.is_initialized())
        throw LIBRPA_INVALID_ARGUMENT(
            label + ": distributed contraction requires an initialized communicator");

    std::string local_error;
    std::size_t npair_size = 0;
    std::size_t result_size = 0;
    if (n_orbitals <= 0 || D_by_atom.empty())
        local_error = label + ": distributed contraction received an empty vertex";
    else if (D_by_atom.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        local_error = label + ": D support is too large for MPI ownership checks";
    else
    {
        npair_size = static_cast<std::size_t>(n_orbitals) * static_cast<std::size_t>(n_orbitals);
        if (npair_size > static_cast<std::size_t>(std::numeric_limits<int>::max()))
            local_error = label + ": ordered-pair dimension exceeds MPI capacity";
        else if (npair_size > std::numeric_limits<std::size_t>::max() / npair_size)
            local_error = label + ": result dimension overflows size_t";
        else
        {
            result_size = npair_size * npair_size;
            if (result_size > static_cast<std::size_t>(std::numeric_limits<int>::max()))
                local_error = label + ": result is too large for MPI_Allreduce";
        }
    }

    // Synchronize the scalar parameter before constructing a result whose
    // shape is used by the next collective.  A malformed rank therefore
    // cannot enter a differently sized allreduce.
    int min_norbitals = 0;
    int max_norbitals = 0;
    MPI_Allreduce(&n_orbitals, &min_norbitals, 1, MPI_INT, MPI_MIN, comm_h.comm);
    MPI_Allreduce(&n_orbitals, &max_norbitals, 1, MPI_INT, MPI_MAX, comm_h.comm);
    if (min_norbitals != max_norbitals)
        local_error = label + ": n_orbitals differs between MPI ranks";

    if (local_error.empty())
    {
        const int npair = static_cast<int>(npair_size);
        for (const auto &[I, D] : D_by_atom)
        {
            if (I < 0 || D.nr != npair || D.nc <= 0)
            {
                local_error = label + ": invalid D atom block shape";
                break;
            }
            try
            {
                validate_finite_matrix(D, label + " D block");
            }
            catch (const std::exception &error)
            {
                local_error = error.what();
                break;
            }
        }
    }
    if (local_error.empty())
    {
        for (const auto &[I, J_W] : W_local_upper)
        {
            if (I < 0)
            {
                local_error = label + ": W atom I is negative";
                break;
            }
            for (const auto &[J, W] : J_W)
            {
                if (J < 0 || I > J)
                {
                    local_error =
                        label + ": W local blocks must use upper-triangular I<=J orientation";
                    break;
                }
                if (W.nr <= 0 || W.nc <= 0)
                {
                    local_error = label + ": W block has an invalid shape";
                    break;
                }
                try
                {
                    validate_finite_matrix(W, label + " W block");
                }
                catch (const std::exception &error)
                {
                    local_error = error.what();
                    break;
                }
                const auto it_I = D_by_atom.find(I);
                const auto it_J = D_by_atom.find(J);
                if (it_I != D_by_atom.end() && it_J != D_by_atom.end() &&
                    (W.nr != it_I->second.nc || W.nc != it_J->second.nc))
                {
                    local_error = label + ": W dimensions do not match D blocks";
                    break;
                }
                if (I == J)
                {
                    if (W.nr != W.nc)
                    {
                        local_error = label + ": diagonal W block is not square";
                        break;
                    }
                    double scale = 1.0;
                    double residual = 0.0;
                    for (int row = 0; row != W.nr; ++row)
                        for (int col = 0; col != W.nc; ++col)
                        {
                            scale = std::max(scale, std::abs(W(row, col)));
                            residual =
                                std::max(residual, std::abs(W(row, col) - std::conj(W(col, row))));
                        }
                    if (!std::isfinite(residual) || residual > 1.0e-10 * scale)
                    {
                        local_error = label + ": diagonal W block is not Hermitian";
                        break;
                    }
                }
            }
            if (!local_error.empty()) break;
        }
    }
    collective_stage_error(comm_h, label + " initial validation", local_error);

    // D support is replicated by the vertex.  Verify that assumption before
    // forming the required-pair list; otherwise different ranks could issue
    // MPI collectives with different vector lengths.
    std::vector<int> local_atoms;
    local_atoms.reserve(D_by_atom.size());
    for (const auto &[I, D] : D_by_atom)
    {
        (void)D;
        local_atoms.push_back(I);
    }
    std::sort(local_atoms.begin(), local_atoms.end());
    const int local_atom_count = static_cast<int>(local_atoms.size());
    std::vector<int> atom_counts(static_cast<std::size_t>(comm_h.nprocs), 0);
    MPI_Allgather(&local_atom_count, 1, MPI_INT, atom_counts.data(), 1, MPI_INT, comm_h.comm);
    std::vector<int> atom_displacements(atom_counts.size(), 0);
    int total_atoms = 0;
    bool atom_count_overflow = false;
    for (std::size_t rank = 0; rank != atom_counts.size(); ++rank)
    {
        if (atom_counts[rank] < 0 ||
            total_atoms > std::numeric_limits<int>::max() - atom_counts[rank])
        {
            atom_count_overflow = true;
            break;
        }
        atom_displacements[rank] = total_atoms;
        total_atoms += atom_counts[rank];
    }
    local_error.clear();
    if (atom_count_overflow) local_error = label + ": D support metadata exceeds MPI capacity";
    collective_stage_error(comm_h, label + " support-count validation", local_error);

    std::vector<int> all_atoms(static_cast<std::size_t>(total_atoms), 0);
    int mpi_dummy_atom = 0;
    const int *local_atoms_buffer = local_atoms.empty() ? &mpi_dummy_atom : local_atoms.data();
    int *all_atoms_buffer = all_atoms.empty() ? &mpi_dummy_atom : all_atoms.data();
    MPI_Allgatherv(local_atoms_buffer, local_atom_count, MPI_INT, all_atoms_buffer,
                   atom_counts.data(), atom_displacements.data(), MPI_INT, comm_h.comm);
    std::set<int> union_atoms(all_atoms.begin(), all_atoms.end());
    local_error.clear();
    for (const int atom : union_atoms)
    {
        if (D_by_atom.count(atom) == 0)
        {
            local_error = label + ": D atom support differs between MPI ranks";
            break;
        }
    }
    collective_stage_error(comm_h, label + " support consensus", local_error);

    std::vector<int> support_atoms(union_atoms.begin(), union_atoms.end());
    const std::size_t support_size = support_atoms.size();
    local_error.clear();
    if (support_size > 0 &&
        support_size > (std::numeric_limits<std::size_t>::max() - support_size) / 2)
        local_error = label + ": required-pair count overflows size_t";
    const std::size_t required_pair_count =
        local_error.empty() ? support_size * (support_size + 1) / 2 : 0;
    if (local_error.empty() &&
        required_pair_count > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        local_error = label + ": too many atom pairs for MPI ownership checks";
    collective_stage_error(comm_h, label + " required-pair capacity", local_error);

    std::vector<std::pair<int, int>> required_pairs;
    required_pairs.reserve(required_pair_count);
    for (std::size_t i = 0; i != support_atoms.size(); ++i)
        for (std::size_t j = i; j != support_atoms.size(); ++j)
            required_pairs.emplace_back(support_atoms[i], support_atoms[j]);

    ComplexMatrix local_result(static_cast<int>(npair_size), static_cast<int>(npair_size));
    local_error.clear();
    if (!W_local_upper.empty())
    {
        try
        {
            local_result = contract_ordered_pair_vertex(D_by_atom, W_local_upper, n_orbitals);
        }
        catch (const std::exception &error)
        {
            local_error = error.what();
        }
        catch (...)
        {
            local_error = "unknown local W contraction failure";
        }
    }
    collective_stage_error(comm_h, label + " local contraction", local_error);

    std::vector<int> local_owners(required_pairs.size(), 0);
    std::vector<int> global_owners(required_pairs.size(), 0);
    for (std::size_t ipair = 0; ipair != required_pairs.size(); ++ipair)
    {
        const auto [I, J] = required_pairs[ipair];
        const auto it_I = W_local_upper.find(I);
        if (it_I != W_local_upper.end() && it_I->second.count(J) != 0) local_owners[ipair] = 1;
    }
    MPI_Allreduce(local_owners.data(), global_owners.data(),
                  static_cast<int>(required_pairs.size()), MPI_INT, MPI_SUM, comm_h.comm);
    local_error.clear();
    for (std::size_t ipair = 0; ipair != required_pairs.size(); ++ipair)
    {
        if (global_owners[ipair] != 1)
        {
            const auto [I, J] = required_pairs[ipair];
            std::ostringstream oss;
            oss << label << ": W block (" << I << ',' << J << ") has " << global_owners[ipair]
                << " MPI owners; exactly one is required";
            local_error = oss.str();
            break;
        }
    }
    collective_stage_error(comm_h, label + " ownership validation", local_error);

    (void)result_size;
    MPI_Allreduce(MPI_IN_PLACE, local_result.c, local_result.size, MPI_CXX_DOUBLE_COMPLEX, MPI_SUM,
                  comm_h.comm);
    return local_result;
}

}  // namespace crpa_test
