#pragma once

// MPI related
#ifdef LIBRPA_USE_MPI
#include <mpi.h>
#else

#ifndef MPI_VERSION

// All wrapped in ifndef, considering that the user might define in there own case.
// In that case, it is assumed that the user-defined macros are consistent with those in LibRPA
// Usually it is the case, since everyone shares the purpose to emulate MPI standard while keeping
// their code compilable with serial compiler.
// I do not know better solution at present, so finger crossed when it is not the case.

#ifndef MPI_Comm
#define MPI_Comm int
#endif

#ifndef MPI_Datatype
#define MPI_Datatype int
#endif

#ifndef MPI_INT
#define MPI_INT 4
#endif

#ifndef MPI_FLOAT
#define MPI_FLOAT 4
#endif

#ifndef MPI_DOUBLE
#define MPI_DOUBLE 8
#endif

#ifndef MPI_LONG
#define MPI_LONG 8
#endif

#ifndef MPI_C_FLOAT_COMPLEX
#define MPI_C_FLOAT_COMPLEX 8
#endif

#ifndef MPI_C_DOUBLE_COMPLEX
#define MPI_C_DOUBLE_COMPLEX 16
#endif

#endif

#endif
