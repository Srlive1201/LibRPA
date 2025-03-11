#pragma once

#include <map>
#include <string>

#include "matrix_m.h"

bool convert_csc(const std::string& filePath, std::map<std::string, Matz>& matrices, std::string& key);
