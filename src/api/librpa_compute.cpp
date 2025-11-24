// Public headers
#include "../../include/librpa_compute.h"

// Internal headers
#include "../core/dataset.h"
#include "../global/instance_manager.h"

// C++ APIs
// The functions here should be just wrappers of C APIs.
namespace librpa
{
    bool test_handler(const LibrpaHandler *h)
    {
        librpa_int::Dataset* p = librpa_int::global::get_dataset_instance(h);
        return is_null_dataset_ptr(p);
    }
}
