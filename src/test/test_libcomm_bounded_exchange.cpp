#include <mpi.h>

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include "Comm/Comm_Trans/Comm_Trans.h"

namespace
{

using Payload = std::vector<double>;
using SendData = std::map<int, std::map<int, Payload>>;
using RecvKey = std::tuple<int, int>;
using RecvData = std::map<RecvKey, Payload>;
using Trans = Comm::Comm_Trans<RecvKey, Payload, SendData, RecvData>;

std::size_t get_payload_scale()
{
    const char *value = std::getenv("LIBCOMM_TEST_DOUBLES");
    if (value == nullptr) return 2000;
    return std::stoull(value);
}

void configure_trans(Trans &trans, const bool copy_merge)
{
    trans.traverse_isend = [](const SendData &data, const int rank_isend,
                              std::function<void(const RecvKey &, const Payload &)> &emit)
    {
        const auto dest_it = data.find(rank_isend);
        if (dest_it == data.end()) return;
        for (const auto &[source_rank, payload] : dest_it->second)
            emit(std::make_tuple(source_rank, rank_isend), payload);
    };
    trans.set_value_recv = [](RecvKey &&key, Payload &&payload, RecvData &data)
    {
        const auto [it, inserted] = data.emplace(std::move(key), std::move(payload));
        if (!inserted) throw std::runtime_error("duplicate received key");
    };
    if (copy_merge)
    {
        trans.flag_lock_set_value = Comm::Comm_Tools::Lock_Type::Copy_merge;
        trans.init_datas_local = [](const int) { return RecvData{}; };
        trans.add_datas = [](RecvData &&source, RecvData &dest)
        {
            for (auto &[key, payload] : source)
            {
                const auto [it, inserted] = dest.emplace(key, std::move(payload));
                if (!inserted) throw std::runtime_error("duplicate merged key");
            }
        };
    }
    else
    {
        trans.flag_lock_set_value = Comm::Comm_Tools::Lock_Type::Lock_free;
    }
}

bool default_uses_bounded_ring(const int rank)
{
    ::unsetenv("LIBCOMM_TRANS_MODE");
    ::setenv("LIBCOMM_TRANS_TRACE", "1", 1);

    SendData send;
    send[rank][rank] = Payload{static_cast<double>(rank)};
    RecvData recv;
    Trans trans(MPI_COMM_WORLD);
    configure_trans(trans, false);

    std::ostringstream captured;
    auto *original = std::cerr.rdbuf(captured.rdbuf());
    trans.communicate(send, recv);
    std::cerr.rdbuf(original);
    ::unsetenv("LIBCOMM_TRANS_TRACE");

    return recv.size() == 1 && captured.str().find("LIBCOMM_TRANS_RING") != std::string::npos;
}

bool rejects_unknown_mode()
{
    ::setenv("LIBCOMM_TRANS_MODE", "not-a-mode", 1);
    SendData send;
    RecvData recv;
    Trans trans(MPI_COMM_WORLD);
    configure_trans(trans, false);
    try
    {
        trans.communicate(send, recv);
    }
    catch (const std::invalid_argument &)
    {
        return true;
    }
    return false;
}

bool rejects_inconsistent_mode(const int rank, const int size)
{
    if (size < 2) return true;
    ::setenv("LIBCOMM_TRANS_MODE", rank == 0 ? "nonblocking" : "sendrecv_ring", 1);
    SendData send;
    RecvData recv;
    Trans trans(MPI_COMM_WORLD);
    configure_trans(trans, false);
    try
    {
        trans.communicate(send, recv);
    }
    catch (const std::invalid_argument &)
    {
        return true;
    }
    return false;
}

bool rejects_inconsistent_chunk_size(const int rank, const int size)
{
    if (size < 2) return true;
    ::setenv("LIBCOMM_TRANS_MODE", "sendrecv_ring", 1);
    ::setenv("LIBCOMM_TRANS_CHUNK_BYTES", rank == 0 ? "1024" : "2048", 1);
    SendData send;
    RecvData recv;
    Trans trans(MPI_COMM_WORLD);
    configure_trans(trans, false);
    try
    {
        trans.communicate(send, recv);
    }
    catch (const std::invalid_argument &)
    {
        return true;
    }
    return false;
}

bool rejects_invalid_chunk_size(const char *value)
{
    ::setenv("LIBCOMM_TRANS_MODE", "sendrecv_ring", 1);
    ::setenv("LIBCOMM_TRANS_CHUNK_BYTES", value, 1);
    SendData send;
    RecvData recv;
    Trans trans(MPI_COMM_WORLD);
    configure_trans(trans, false);
    try
    {
        trans.communicate(send, recv);
    }
    catch (const std::invalid_argument &)
    {
        return true;
    }
    return false;
}

bool checks_empty_ring_exchange()
{
    ::setenv("LIBCOMM_TRANS_MODE", "sendrecv_ring", 1);
    ::setenv("LIBCOMM_TRANS_CHUNK_BYTES", "1024", 1);
    SendData send;
    RecvData recv;
    Trans trans(MPI_COMM_WORLD);
    configure_trans(trans, false);
    trans.communicate(send, recv);
    return recv.empty();
}

bool checks_single_traversal_per_destination(const int size)
{
    ::setenv("LIBCOMM_TRANS_MODE", "sendrecv_ring", 1);
    ::setenv("LIBCOMM_TRANS_CHUNK_BYTES", "1024", 1);
    SendData send;
    RecvData recv;
    Trans trans(MPI_COMM_WORLD);
    configure_trans(trans, false);
    std::vector<int> traversals(static_cast<std::size_t>(size), 0);
    trans.traverse_isend = [&traversals](
        const SendData &, const int destination,
        std::function<void(const RecvKey &, const Payload &)> &)
    {
        ++traversals.at(static_cast<std::size_t>(destination));
    };
    trans.communicate(send, recv);
    return std::all_of(traversals.begin(), traversals.end(),
                       [](const int count) { return count == 1; });
}

bool rejects_rank_local_serialization_failure(const int rank)
{
    ::setenv("LIBCOMM_TRANS_MODE", "sendrecv_ring", 1);
    ::setenv("LIBCOMM_TRANS_CHUNK_BYTES", "1024", 1);
    SendData send;
    send[rank][rank] = Payload{static_cast<double>(rank)};
    RecvData recv;
    Trans trans(MPI_COMM_WORLD);
    configure_trans(trans, false);
    trans.traverse_isend = [rank](const SendData &, const int,
                                  std::function<void(const RecvKey &, const Payload &)> &)
    {
        if (rank == 0) throw std::runtime_error("injected serialization failure");
    };
    try
    {
        trans.communicate(send, recv);
    }
    catch (const std::runtime_error &)
    {
        return true;
    }
    return false;
}

bool rejects_rank_local_deserialization_failure(const int rank)
{
    ::setenv("LIBCOMM_TRANS_MODE", "sendrecv_ring", 1);
    ::setenv("LIBCOMM_TRANS_CHUNK_BYTES", "1024", 1);
    SendData send;
    send[rank][rank] = Payload{static_cast<double>(rank)};
    RecvData recv;
    Trans trans(MPI_COMM_WORLD);
    configure_trans(trans, false);
    trans.set_value_recv = [rank](RecvKey &&key, Payload &&payload, RecvData &data)
    {
        if (rank == 0) throw std::runtime_error("injected deserialization failure");
        data.emplace(std::move(key), std::move(payload));
    };
    try
    {
        trans.communicate(send, recv);
    }
    catch (const std::runtime_error &)
    {
        return true;
    }
    return false;
}

int check_exchange(const char *mode, const bool copy_merge, const int rank, const int size)
{
    ::setenv("LIBCOMM_TRANS_MODE", mode, 1);
    if (std::string(mode) == "sendrecv_ring")
        ::setenv("LIBCOMM_TRANS_CHUNK_BYTES", "1024", 1);
    else
        ::unsetenv("LIBCOMM_TRANS_CHUNK_BYTES");
    const std::size_t scale = get_payload_scale();

    SendData send;
    for (int dest = 0; dest != size; ++dest)
    {
        const std::size_t count =
            scale * static_cast<std::size_t>(rank + 1) * static_cast<std::size_t>(dest + 1);
        Payload payload(count);
        for (std::size_t i = 0; i != count; ++i)
            payload[i] = 1000000.0 * rank + 1000.0 * dest + static_cast<double>(i);
        send[dest][rank] = std::move(payload);
    }

    RecvData recv;
    Trans trans(MPI_COMM_WORLD);
    configure_trans(trans, copy_merge);
    trans.communicate(send, recv);

    int failures = recv.size() == static_cast<std::size_t>(size) ? 0 : 1;
    for (int source = 0; source != size; ++source)
    {
        const RecvKey key{source, rank};
        const auto it = recv.find(key);
        if (it == recv.end())
        {
            ++failures;
            continue;
        }
        const std::size_t expected_count =
            scale * static_cast<std::size_t>(source + 1) * static_cast<std::size_t>(rank + 1);
        if (it->second.size() != expected_count)
        {
            ++failures;
            continue;
        }
        for (std::size_t i : {std::size_t{0}, expected_count / 2, expected_count - 1})
        {
            const double expected = 1000000.0 * source + 1000.0 * rank + static_cast<double>(i);
            if (it->second[i] != expected) ++failures;
        }
    }
    return failures;
}

}  // namespace

int main(int argc, char *argv[])
{
    int provided = 0;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    if (provided < MPI_THREAD_FUNNELED)
    {
        MPI_Finalize();
        return 1;
    }

    int rank = 0;
    int size = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    const char *backend_env = std::getenv("LIBCOMM_TEST_BACKEND");
    const std::string backend = backend_env == nullptr ? "sendrecv_ring" : backend_env;
    int local_failures = default_uses_bounded_ring(rank) ? 0 : 1;
    local_failures += rejects_unknown_mode() ? 0 : 1;
    local_failures += rejects_inconsistent_mode(rank, size) ? 0 : 1;
    local_failures += rejects_inconsistent_chunk_size(rank, size) ? 0 : 1;
    local_failures += rejects_invalid_chunk_size("0") ? 0 : 1;
    local_failures += rejects_invalid_chunk_size("not-a-number") ? 0 : 1;
    local_failures += checks_empty_ring_exchange() ? 0 : 1;
    local_failures += checks_single_traversal_per_destination(size) ? 0 : 1;
    local_failures += rejects_rank_local_serialization_failure(rank) ? 0 : 1;
    local_failures += rejects_rank_local_deserialization_failure(rank) ? 0 : 1;
    if (backend == "nonblocking" || backend == "all")
    {
        local_failures += check_exchange("nonblocking", false, rank, size);
        local_failures += check_exchange("nonblocking", true, rank, size);
    }
    if (backend == "sendrecv_ring" || backend == "all")
    {
        local_failures += check_exchange("sendrecv_ring", false, rank, size);
        local_failures += check_exchange("sendrecv_ring", true, rank, size);
    }
    if (backend != "nonblocking" && backend != "sendrecv_ring" && backend != "all")
        ++local_failures;

    int global_failures = 0;
    MPI_Allreduce(&local_failures, &global_failures, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0)
        std::cout << "LIBCOMM_BOUNDED_EXCHANGE backend=" << backend << " size=" << size
                  << " failures=" << global_failures << std::endl;

    MPI_Finalize();
    return global_failures == 0 ? 0 : 2;
}
