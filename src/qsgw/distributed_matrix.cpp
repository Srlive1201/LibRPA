#include "distributed_matrix.h"

#include "../math/scalapack_connector.h"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace librpa_int
{
namespace qsgw
{

Matz collect_blacs_matrix_root(const Matz& local,
                               const ArrayDesc& distributed_descriptor)
{
    if (!distributed_descriptor.is_initialized())
    {
        throw std::invalid_argument(
            "QSGW distributed matrix descriptor is not initialized");
    }
    if (distributed_descriptor.m() <= 0 ||
        distributed_descriptor.n() <= 0)
    {
        throw std::invalid_argument(
            "QSGW distributed matrix must have positive global dimensions");
    }
    if (local.major() != MAJOR::COL)
    {
        throw std::invalid_argument(
            "QSGW distributed matrix must use column-major local storage");
    }
    if (local.nr() != distributed_descriptor.m_loc() ||
        local.nc() != distributed_descriptor.n_loc())
    {
        throw std::invalid_argument(
            "QSGW local matrix shape does not match its BLACS descriptor");
    }

    ArrayDesc root_descriptor(distributed_descriptor.ictxt());
    root_descriptor.init(distributed_descriptor.m(),
                         distributed_descriptor.n(),
                         distributed_descriptor.m(),
                         distributed_descriptor.n(),
                         distributed_descriptor.irsrc(),
                         distributed_descriptor.icsrc());

    Matz source_dummy(1, 1, MAJOR::COL);
    const cplxdb* source = local.nr() > 0 && local.nc() > 0
                                ? local.ptr()
                                : source_dummy.ptr();
    Matz transfer_buffer = root_descriptor.is_src()
                               ? Matz(distributed_descriptor.m(),
                                      distributed_descriptor.n(), MAJOR::COL)
                               : Matz(1, 1, MAJOR::COL);
    ScalapackConnector::pgemr2d_f(
        distributed_descriptor.m(), distributed_descriptor.n(),
        source, 1, 1, distributed_descriptor.desc,
        transfer_buffer.ptr(), 1, 1, root_descriptor.desc,
        distributed_descriptor.ictxt());

    if (root_descriptor.is_src())
    {
        return transfer_buffer;
    }
    return Matz(0, 0, MAJOR::COL);
}

void broadcast_spin_k_matrix_map(SpinKMatrixMap& values,
                                 const int root,
                                 const MpiCommHandler& communicator)
{
    communicator.check_initialized();
    if (root < 0 || root >= communicator.nprocs)
    {
        throw std::invalid_argument("QSGW matrix broadcast root is invalid");
    }

    int valid = 1;
    int block_count = 0;
    std::vector<int> metadata;
    if (communicator.myid == root)
    {
        if (values.empty()) valid = 0;
        for (const auto& [spin, by_kpoint] : values)
        {
            if (by_kpoint.empty()) valid = 0;
            for (const auto& [kpoint, matrix] : by_kpoint)
            {
                if (matrix.nr() <= 0 || matrix.nc() <= 0 ||
                    matrix.nr() > std::numeric_limits<int>::max() /
                                      matrix.nc())
                {
                    valid = 0;
                    continue;
                }
                for (int index = 0; index < matrix.nr() * matrix.nc(); ++index)
                {
                    const cplxdb value = matrix.ptr()[index];
                    if (!std::isfinite(value.real()) ||
                        !std::isfinite(value.imag()))
                    {
                        valid = 0;
                    }
                }
                metadata.insert(metadata.end(), {
                    spin, kpoint, matrix.nr(), matrix.nc(),
                    static_cast<int>(matrix.major())});
                ++block_count;
            }
        }
        if (block_count <= 0) valid = 0;
    }

    communicator.bcast(&valid, 1, root);
    if (!valid)
    {
        throw std::invalid_argument(
            "QSGW root matrix map is empty, malformed, or non-finite");
    }
    communicator.bcast(&block_count, 1, root);
    if (communicator.myid != root)
    {
        metadata.resize(static_cast<std::size_t>(5 * block_count));
    }
    communicator.bcast(metadata.data(), 5 * block_count, root);

    SpinKMatrixMap received;
    for (int block = 0; block < block_count; ++block)
    {
        const int offset = 5 * block;
        const int spin = metadata[offset];
        const int kpoint = metadata[offset + 1];
        const int rows = metadata[offset + 2];
        const int columns = metadata[offset + 3];
        const int major_raw = metadata[offset + 4];
        if (rows <= 0 || columns <= 0 ||
            (major_raw != static_cast<int>(MAJOR::ROW) &&
             major_raw != static_cast<int>(MAJOR::COL)))
        {
            throw std::runtime_error(
                "QSGW matrix broadcast metadata is invalid");
        }
        Matz matrix(rows, columns, static_cast<MAJOR>(major_raw));
        if (communicator.myid == root)
        {
            const auto spin_it = values.find(spin);
            if (spin_it == values.end())
                throw std::runtime_error(
                    "QSGW root matrix map changed during broadcast");
            const auto kpoint_it = spin_it->second.find(kpoint);
            if (kpoint_it == spin_it->second.end())
                throw std::runtime_error(
                    "QSGW root matrix map changed during broadcast");
            matrix = kpoint_it->second.copy();
        }
        communicator.bcast(matrix.ptr(), rows * columns, root);
        if (!received[spin].emplace(kpoint, std::move(matrix)).second)
        {
            throw std::runtime_error(
                "QSGW matrix broadcast contains duplicate keys");
        }
    }
    values = std::move(received);
}

} // namespace qsgw
} // namespace librpa_int
