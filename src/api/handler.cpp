// Public headers (prefixed by librpa)
#include "../../include/librpa_handler.h"

// Internal headers
#include "../global/instance_manager.h"

// C APIs
LibrpaHandler* librpa_create_handler(int comm)
{
    return librpa_int::global::push_back_dataset(comm);
}

void librpa_destroy_handler(LibrpaHandler *h)
{
    if (!h) return;
    librpa_int::global::destroy_dataset(h);
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

