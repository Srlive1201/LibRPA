#pragma once
#include <vector>

// Public headers
#include "../../include/librpa_handler.h"

// Internal headers
#include "../core/dataset.h"

namespace librpa_int
{

extern std::vector<Dataset*> manager;

Dataset* get_dataset_instance(const LibrpaHandler *h);

}
