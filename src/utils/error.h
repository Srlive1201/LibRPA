#pragma once

#include <string>
#include <stdexcept>

namespace librpa_int
{

#define LIBRPA_ERROR_PREFIX std::string(__FILE__) + ":" + \
    std::to_string(__LINE__) + ":" + std::string(__FUNCTION__) + ": "

#define LIBRPA_RUNTIME_ERROR(msg) std::runtime_error(LIBRPA_ERROR_PREFIX + msg)

}
