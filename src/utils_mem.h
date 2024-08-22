#pragma once

namespace LIBRPA
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


} /* end of namespace utils */

} /* end of namespace LIBRPA */

