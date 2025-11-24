#pragma once

namespace librpa_int
{

/*!
 * @brief Create directory and its parent directories if necessary.
 *
 * If the directory exists already, the function returns silently.
 * Need to barrier around this function if processes other than root process
 * wants to write files under the directory after creation.
 *
 * @param[in]  dname         Name of the directory to create
 * @param[in]  root_process  Process filter.
 *                           Only the processes that parses 0 to root_process will create the target directory.
 */
void create_directories(const char *dname, int root_process);

}
