#pragma once
/**
 * @file librpa_global.h
 * @brief Global functions for LibRPA.
 *
 * These functions include:
 * - initialization and finalization of the global LibRPA environment,
 * - version check,
 * - output of profiling statistics,
 * - sanity check test function for debug use (`librpa_test`).
 *
 * Global functions do not require initialization before hand.
 */

#include "librpa_enums.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Get build information string.
 *
 * @return C-string containing build information (e.g., compiler, options, date).
 */
const char* librpa_get_build_info(void);

/**
 * @brief Get major version number.
 *
 * @return Major version component (X in X.Y.Z).
 */
int librpa_get_major_version(void);

/**
 * @brief Get minor version number.
 *
 * @return Minor version component (Y in X.Y.Z).
 */
int librpa_get_minor_version(void);

/**
 * @brief Get patch version number.
 *
 * @return Patch version component (Z in X.Y.Z).
 */
int librpa_get_patch_version(void);

/**
 * @brief Initialize the global LibRPA environment.
 *
 * Must be called after MPI_Init() and before any other LibRPA functions.
 * Sets up logging, profiling, and parallel communication infrastructure.
 *
 * @param switch_redirect_stdout If LIBRPA_SWITCH_ON, redirect stdout to a file.
 *                               Default is LIBRPA_SWITCH_OFF.
 * @param redirect_path          Path for redirected output (default: "stdout").
 *                               Only used when switch_redirect_stdout is ON.
 * @param switch_process_output  If LIBRPA_SWITCH_ON, enable per-process output files.
 *                               Default is LIBRPA_SWITCH_ON.
 */
void librpa_init_global(LibrpaSwitch switch_redirect_stdout = LIBRPA_SWITCH_OFF,
                        const char *redirect_path = "stdout",
                        LibrpaSwitch switch_process_output = LIBRPA_SWITCH_ON);

/**
 * @brief Finalize the global LibRPA environment.
 *
 * Must be called after all LibRPA operations are complete and before MPI_Finalize().
 * Releases all global resources.
 */
void librpa_finalize_global(void);

/**
 * @brief Run internal self-tests.
 *
 * Performs basic sanity checks on LibRPA functionality.
 * Useful for debugging.
 */
void librpa_test(void);

/**
 * @brief Print profiling information.
 *
 * Outputs timing and memory usage statistics collected during computation.
 * Only produces meaningful output if profiling was enabled at compile time.
 */
void librpa_print_profile(void);

#ifdef __cplusplus
}
#endif
