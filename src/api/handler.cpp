// Public headers (prefixed by librpa)
#include "../../include/librpa_handler.h"

// Internal headers
#include "instance_manager.h"

// C APIs
LibrpaHandler* librpa_create_handler(int comm)
{
    return librpa_int::api::push_back_dataset(comm);
}

void librpa_destroy_handler(LibrpaHandler *h)
{
    if (!h) return;
    librpa_int::api::destroy_dataset(h);
    delete h;
}
