#include "instance_manager.h"

#include <cstddef>
#include "librpa.hpp"

namespace librpa_int
{

namespace api
{

// manager[0] is a reserved sentinel (always nullptr); real instances start at 1.
std::vector<Dataset *> manager{nullptr};

librpa_int::Dataset* get_dataset_instance(const LibrpaHandler *h)
{
    librpa_int::Dataset* p = nullptr;
    const int total_size = manager.size();
    if (h != nullptr)
    {
        const auto id = h->instance_id_;
        if (id >= 0 && id < total_size) p = manager[id];
    }
    return p;
}

librpa_int::Dataset* get_dataset_instance(const librpa::Handler &h)
{
    return get_dataset_instance(h.get_c_handler());
}

LibrpaHandler* push_back_dataset(int comm)
{
    // create a new instance and append it to the manager
    int instance_id = manager.size();
    librpa_int::Dataset *obj = new librpa_int::Dataset(comm);
    manager.emplace_back(obj);

    // initialize a binding handler
    LibrpaHandler* h = new LibrpaHandler {instance_id};
    return h;
}

void destroy_dataset(LibrpaHandler* h)
{
    if (!h) return;
    // destroy the instance
    if (h->instance_id_ > 0)
    {
        const auto id = h->instance_id_;
        const int total_size = manager.size();
        // Invalid handler that was manually created with hand-picked id,
        // either oversubscription
        if (id >= total_size) return;
        auto &p = manager[id];
        // or pointed to an already released instance
        if (!p) return;
        // free the instance and point it to null pointer
        delete p;
        manager[id] = nullptr;
    }
}

} /* namespace api */

} /* namespace librpa_int */
