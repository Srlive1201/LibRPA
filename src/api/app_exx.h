/*!
 * @file      app_exx.h
 * @author    Min-Ye Zhang
 * @date      2024-07-02
 */
#pragma once
#include <vector>

namespace LIBRPA
{

namespace app
{

/*!
 * @brief Compute the exact exchange (EXX) energy for states at specified k-points
 *
 * @param[in] i_state_low       The lowest index of state (included) to compute
 * @param[in] i_state_high      The highest index of state (excluded) to compute
 * @param[in] n_kpoints_task    The number of k-points to return in the called process.
 *                              When equal to 0, an empty vector will be returned.
 *                              When less than 0, all k-points will be computed.
 *                              Otherwise, the states at k-points whose indices
 *                              are stored in `i_kpoints_task` will be computed.
 * @param[in] i_kpoints_task    The indices of k-points to compute EXX energy.
 *
 * @returns   `std::vector<double>`, exchange energy of states.
 * If `i_state_low` is no less than `i_state_high`, or `n_kpoints_task` is 0, the vector is empty.
 * Otherwise the vector is of size `n_spins` * `n_kpoints_compute` * (`i_state_high` - `i_state_low`),
 * where `n_kpoints_compute` equals to total number of k-points if `n_kpoints_task` < 1, and
 * `n_kpoints_task` otherwise. The indices runs in the order of states, k-points and spins.
 */
std::vector<double>
compute_exx_orbital_energy_(int i_state_low, int i_state_high,
                            int n_kpoints_task, const int *i_kpoints_task);

} /* end of namespace app */

} /* end of namespace LIBRPA */
