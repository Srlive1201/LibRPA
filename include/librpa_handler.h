#pragma once

/**
 * @file librpa_handler.h
 * @brief Handler management for LibRPA instances.
 */

// C APIs
#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Opaque handle to a LibRPA instance.
 *
 * This structure encapsulates the internal state of a LibRPA computation.
 * Users should not access its members directly; use the provided API functions.
 *
 * @note Do not create this structure manually. Use librpa_create_handler() instead.
 */
typedef struct
{
    /** Internal instance ID used by the manager. Do not modify manually. */
    const int instance_id_;
} LibrpaHandler;

/**
 * @brief Create a new LibRPA handler instance.
 *
 * Allocates and initializes a new LibRPA handler associated with the given
 * MPI communicator. This handle is used in subsequent API calls to
 * perform RPA/EXX/G0W0 calculations.
 *
 * @param comm MPI communicator (e.g., MPI_COMM_WORLD).
 *
 * @return Pointer to newly created LibrpaHandler. Must be freed with
 *         librpa_destroy_handler() when no longer needed.
 *
 * @see librpa_destroy_handler
 */
LibrpaHandler* librpa_create_handler(int comm);

/**
 * @brief Destroy a LibRPA handler instance.
 *
 * Frees all internal resources associated with the handler. After this call,
 * the handler pointer becomes invalid and should not be used.
 *
 * @param h Pointer to the LibrpaHandler to destroy. If NULL, no action is taken.
 *
 * @see librpa_create_handler
 */
void librpa_destroy_handler(LibrpaHandler *h);

#ifdef __cplusplus
}
#endif
