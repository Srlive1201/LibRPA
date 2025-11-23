#pragma once

// MPI related
#ifdef __USE_MPI__
#include <mpi.h>
#else

// All wrapped in ifndef, considering that the user might define in there own case.
// In that case, it is assumed that the user-defined macros are consistent with those in LibRPA
// Usually it is the case, since everyone shares the purpose to emulate MPI standard while keeping
// their code compilable with serial compiler.
// I do not know better solution at present, so finger crossed when it is not the case.

#ifndef MPI_Comm
#define MPI_Comm int
#endif

#endif
