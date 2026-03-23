#pragma once
/**
 * @file librpa.h
 * @brief Main include file for LibRPA.
 *
 * This is the primary header file for LibRPA. It aggregates all public APIs:
 * - Global environment functions (librpa_global.h)
 * - Runtime options (librpa_options.h)
 * - Handler management (librpa_handler.h)
 * - Input data APIs (librpa_input.h)
 * - Computation APIs (librpa_compute.h)
 * - C++ wrapper (librpa.hpp)
 *
 * Usage:
 * @code
 * #include "librpa.h"
 * @endcode
 */

#include "librpa_global.h"

// Input data APIs
#include "librpa_input.h"

// Computation APIs
#include "librpa_compute.h"

// C++ APIs: wrappers to the stable C APIs
#ifdef __cplusplus
#include "librpa.hpp"
#endif
