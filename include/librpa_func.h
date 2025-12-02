#pragma once
// Macros for function declaration and definition

// C functions for handle objects with and without parameters, without runtime options
#define LIBRPA_C_H_FUNC_WRAP(ret, func, ...) ret func(LibrpaHandler *h, __VA_ARGS__)
#define LIBRPA_C_H_FUNC_WRAP_NOPAR(ret, func) ret func(LibrpaHandler *h)

// C functions for handle objects with and without parameters, with runtime options
#define LIBRPA_C_H_FUNC_WRAP_WOPT(ret, func, ...) ret func(LibrpaHandler *h, const LibrpaOptions *opts, __VA_ARGS__)
#define LIBRPA_C_H_FUNC_WRAP_WOPT_NOPAR(ret, func) ret func(LibrpaHandler *h, const LibrpaOptions *opts)
