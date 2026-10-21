#pragma once
#include "../core/qsgw/task.h"
#include "librpa_qsgw.h"
namespace librpa_int::api
{
librpa_int::qsgw::Task& qsgw_task_instance(LibrpaQsgwTask*);
}
