#pragma once

#include <exception>
#include <string>

#include "../io/global_io.h"
#include "../utils/error.h"
#include "base_mpi.h"

namespace librpa_int
{

// Propagate rank-local preparation errors before entering the next collective.
// The supplied operation must itself preserve collective participation.
template <typename Function>
void collective_all_stage(const MpiCommHandler &comm, const std::string &label,
                          Function &&operation)
{
    std::exception_ptr local_error;
    int failed = 0, any_failed = 0;
    try
    {
        operation();
    }
    catch (const std::exception &error)
    {
        failed = 1;
        local_error = std::current_exception();
        global::lib_printf(LIBRPA_VERBOSE_CRITICAL, "%s failed: %s\n", label.c_str(), error.what());
    }
    catch (...)
    {
        failed = 1;
        local_error = std::current_exception();
    }
    comm.allreduce(&failed, &any_failed, 1, MPI_MAX);
    if (any_failed)
    {
        if (local_error) std::rethrow_exception(local_error);
        throw LIBRPA_RUNTIME_ERROR(label + " failed on another rank");
    }
}

}  // namespace librpa_int
