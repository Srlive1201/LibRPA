#pragma once

#include <ostream>

namespace librpa_int
{

namespace utils
{

/*!
 * @brief Display the free memory of the system.
 */
void display_free_mem();


/*!
 * @brief Release system free memory.
 */
void release_free_mem();


/*!
 * @brief Get the number of total memory on a node
 *
 * @param[out] Total memory available on the node in GB
 *
 * @retval integer, 0 for successful get
 *
 * @note
 *   Currently only works with Linux and BSD operating systems.
 */
int get_node_total_mem(double &total_mem_gb);


/*!
 * @brief Get the number of free memory on a node
 *
 * @param[out] Free memory available on the node in GB
 *
 * @retval integer, 0 for successful get
 *
 * @note
 *   Currently only works with Linux operating system.
 */
int get_node_free_mem(double &free_mem_gb);


/*!
 * @brief Report virtual page usage per process
 *
 * @param os ostream handler that the output should be directed,
 *           favorably each process have its own .
 *
 * @note
 *   Currently only works with Linux operating system.
 */
void report_virtual_pages(std::ostream &os);


} /* end of namespace utils */

} /* end of namespace librpa_int */

