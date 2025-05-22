#ifndef TASK_QSGW_band_H
#define TASK_QSGW_band_H
#include <map>

#include "complexmatrix.h"
#include "vector3_order.h"
// 声明 QSGW 计算任务的函数
void task_qsgw_band(std::map<Vector3_Order<double>, ComplexMatrix>& sinvS);
void plot_homo_lumo_vs_iterations();
#endif  // TASK_QSGW_H
