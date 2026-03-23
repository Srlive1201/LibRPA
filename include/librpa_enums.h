#pragma once

/**
 * @file librpa_enums.h
 * @brief Enumeration types and constants for LibRPA.
 *
 * This file defines the enumeration types used throughout LibRPA for
 * parallel routing, time/frequency grids, and runtime controls.
 */

// C enums
#ifdef __cplusplus
extern "C" {
#endif

/** @brief Undefined or unset value for integer parameters. */
#define LIBRPA_UNSET -101

/** @brief Automatic selection value. LibRPA will choose appropriate setting. */
#define LIBRPA_AUTO -51

/** @brief Switch value for disabled/off state (equivalent to false). */
#define LIBRPA_SWITCH_OFF 0

/** @brief Switch value for enabled/on state (equivalent to true). */
#define LIBRPA_SWITCH_ON 1

/** @brief Verbose level for debug output. */
#define LIBRPA_VERBOSE_DEBUG 4

/** @brief Verbose level for warning messages. */
#define LIBRPA_VERBOSE_WARN 3

/** @brief Verbose level for informational messages. */
#define LIBRPA_VERBOSE_INFO 2

/** @brief Verbose level for critical/error messages only. */
#define LIBRPA_VERBOSE_CRITICAL 1

/** @brief Silent mode - no output. */
#define LIBRPA_VERBOSE_SILENT 0

// Reserved for future DFT code interfaces
// #define LIBRPA_KIND_AIMS 100
// #define LIBRPA_KIND_ABACUS 101
// #define LIBRPA_KIND_OPENMX 102
// #define LIBRPA_KIND_PYSCF 103

/** Number of parallel routing types available. */
#define LIBRPA_ROUTING_COUNT 5

/**
 * @brief Parallel routing strategy for distributed memory calculations.
 *
 * Specifies how the computation is distributed across MPI processes.
 */
typedef enum
{
    ROUTING_UNSET = LIBRPA_UNSET,  ///< Use default routing (unset)
    AUTO = LIBRPA_AUTO,             ///< Automatically select optimal routing
    RTAU = 0,                      ///< Real-space tau (time) decomposition
    ATOMPAIR = 1,                  ///< Atom-pair parallelization
    LIBRI = 2,                     ///< Use LibRI for RI basis operations
} LibrpaParallelRouting;

/** Number of time/frequency grid types available. */
#define LIBRPA_TFGRID_COUNT 7

/**
 * @brief Type of time or frequency grid for integration.
 *
 * Different grid types offer different convergence properties
 * for RPA and GW calculations.
 */
typedef enum
{
    TFGRID_UNSET = LIBRPA_UNSET,   ///< Use default grid (unset)
    GaussLegendre = 0,              ///< Gauss-Legendre quadrature
    GaussChebyshevI = 1,           ///< Gauss-Chebyshev type I
    GaussChebyshevII = 2,          ///< Gauss-Chebyshev type II
    Minimax = 3,                   ///< Minimax grid
    EvenSpaced = 4,                ///< Evenly spaced grid
    EvenSpaced_TF = 5,             ///< Evenly spaced in time-frequency
} LibrpaTimeFreqGrid;

/**
 * @brief Boolean switch type.
 *
 * Use LIBRPA_SWITCH_ON (1) or LIBRPA_SWITCH_OFF (0) to set.
 */
typedef int LibrpaSwitch;

/**
 * @brief Type of DFT code (reserved for future use).
 *
 * Currently reserved for identifying the source of wavefunction data.
 * Not yet implemented.
 */
typedef int LibrpaKind;

/**
 * @brief Verbosity level for runtime output.
 *
 * Controls the amount of information printed during computation.
 * Use one of LIBRPA_VERBOSE_* constants.
 */
typedef int LibrpaVerbose;

#ifdef __cplusplus
}
#endif
