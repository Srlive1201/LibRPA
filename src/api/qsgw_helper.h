#pragma once
#include "dataset.h"
#include "../core/qsgw/task.h"
#include "librpa_qsgw.h"
namespace librpa_int::api
{
librpa_int::qsgw::Task& qsgw_task_instance(LibrpaQsgwTask*);
librpa_int::dataset_ptr_t qsgw_task_dataset(LibrpaQsgwTask*);
}  // namespace librpa_int::api
