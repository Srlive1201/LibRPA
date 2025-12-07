#pragma once
#include <iterator>
#include <vector>
#include <memory>

// Public headers
#include "librpa_handler.h"

// Internal headers
#include "dataset.h"

namespace librpa_int
{

namespace api
{

//! Manager of dataset instances created by user
extern std::vector<dataset_ptr_t> manager;

dataset_ptr_t get_dataset_instance(const LibrpaHandler *h);

LibrpaHandler* push_back_dataset(int comm);

//! Free the data instance that the handler binds
void destroy_dataset(LibrpaHandler* h);

}

}
