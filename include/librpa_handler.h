#pragma once

// C APIs
#ifdef __cplusplus
extern "C" {
#endif

// Do not create by hand
typedef struct
{
    // The only member: an integer that maps to the index of working instance in the manager.
    // The program may not work properly with manually created handler
    const int instance_id_;
} LibrpaHandler;

LibrpaHandler* librpa_create_handler(int comm);

void librpa_destroy_handler(LibrpaHandler *h);

#ifdef __cplusplus
}
#endif

// C++ APIs
#ifdef __cplusplus

namespace librpa
{

class Handler
{
private:
    LibrpaHandler *h;
public:
    Handler(int comm);
    LibrpaHandler *get_c_handler() const { return h; }
    ~Handler();
};

} /* namespace librpa */

#endif
