#ifndef TASK_QSGW_H
#define TASK_QSGW_H

#include <vector>

// 声明全局变量（仅声明，不定义）
extern std::vector<double> homo_values;
extern std::vector<double> lumo_values;
extern std::vector<double> efermi_values;
extern std::vector<int> iteration_numbers;
// 声明 QSGW 计算任务的函数
void task_qsgw();
void plot_homo_lumo_vs_iterations();
#endif // TASK_QSGW_H
