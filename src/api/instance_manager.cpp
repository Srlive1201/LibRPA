#include "instance_manager.h"

#include <cstddef>
#include <memory>
#include "librpa.hpp"

namespace librpa_int
{

namespace api
{

// manager[0] is a reserved sentinel (always nullptr); real instances start at 1.
std::vector<std::shared_ptr<Dataset>> manager{nullptr};

std::shared_ptr<Dataset> get_dataset_instance(const LibrpaHandler *h)
{
    if (h == nullptr) return nullptr;
    return manager.at(h->instance_id_);
}

std::shared_ptr<Dataset> get_dataset_instance(const librpa::Handler &h)
{
    return get_dataset_instance(h.get_c_handler());
}

LibrpaHandler* push_back_dataset(int comm)
{
    // create a new instance and append it to the manager
    int instance_id = manager.size();
    auto p = std::make_shared<Dataset>(comm);
    manager.emplace_back(p);

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
        auto p = manager[id];
        // or pointed to an already released instance
        p.reset();
        // free the instance and point it to null pointer
        manager[id] = nullptr;
    }
}

} /* namespace api */

} /* namespace librpa_int */
