#include "instance_manager.h"

#include <cstddef>

namespace librpa_int
{

// manager[0] is a reserved sentinel (always nullptr); real instances start at 1.
std::vector<Dataset *> manager{nullptr};

librpa_int::Dataset* get_dataset_instance(const LibrpaHandler *h)
{
    librpa_int::Dataset* p = nullptr;
    const int total_size = librpa_int::manager.size();
    if (h != nullptr)
    {
        const auto id = h->instance_id_;
        if (id >= 0 && id < total_size) p = librpa_int::manager[id];
    }
    return p;
}

}
