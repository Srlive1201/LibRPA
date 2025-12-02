#pragma once

#include <string>
#include <sstream>

#include "envs_cmake.h"

namespace librpa_int
{

/*!
 * @brief Print the CMake build information of LibRPA
 */
void print_cmake_info();

/*!
 * @brief Format and return the CMake build information of LibRPA
 */
const std::string& cmake_info_storage();

} /* end of namespace librpa_int */
