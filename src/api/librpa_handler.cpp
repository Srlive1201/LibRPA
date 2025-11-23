// Public headers (prefixed by librpa)
#include "../../include/librpa_handler.h"

// Internal headers
#include "../global/instance_manager.h"

#include "../core/dataset.h"

// C APIs
LibrpaHandler* librpa_create_handler(int comm)
{
    // create a new instance and append it to the manager
    int instance_id = librpa_int::manager.size();
    librpa_int::Dataset *obj = new librpa_int::Dataset(comm);
    librpa_int::manager.emplace_back(obj);

    // initialize a binding handler
    LibrpaHandler* h = new LibrpaHandler {instance_id};
    return h;
}

// free the data instance that the handler binds
static void free_handler_data(LibrpaHandler *h)
{
    if (!h) return;
    // destroy the instance
    if (h->instance_id_ > 0)
    {
        const auto id = h->instance_id_;
        const int total_size = librpa_int::manager.size();
        // Invalid handler that was manually created with hand-picked id,
        // either oversubscription
        if (id >= total_size) return;
        auto &p = librpa_int::manager[id];
        // or pointed to an already released instance
        if (!p) return;
        // free the instance and point it to null pointer
        delete p;
        librpa_int::manager[id] = nullptr;
    }
}

void librpa_destroy_handler(LibrpaHandler *h)
{
    if (!h) return;
    free_handler_data(h);
    delete h;
}

// C++ APIs

namespace librpa
{

Handler::Handler(int comm): h(nullptr)
{
    h = ::librpa_create_handler(comm);
}

Handler::~Handler()
{
    ::librpa_destroy_handler(h);
}

}

