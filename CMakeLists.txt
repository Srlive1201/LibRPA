cmake_minimum_required(VERSION 3.10)

project(LibRPA)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_C_COMPILER "mpiicc" CACHE STRING "" FORCE)
set(CMAKE_CXX_COMPILER "mpiicpc" CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS "-g -O2 -qopenmp" CACHE STRING "" FORCE)

set(LAPACK_DIR "$ENV{MKLROOT}/lib/intel64")
set(LIBS "-lmkl_intel_lp64 -lmkl_sequential -lmkl_core")

set(LAPACK_LIBRARIES "-L${LAPACK_DIR} ${LIBS}")


# generate_library_targets(LIB_PATHS LIBS)

add_executable(chi0_main.exe
        lib_main.cpp
        cal_periodic_chi0.cpp
        aperiodic_chi0.cpp
        Gauss_Quadrature.cpp
        vector3_order.cpp
        matrix.cpp
        complexmatrix.cpp
        matrix3.cpp
        input.cpp
        parallel_mpi.cpp
        profiler.cpp
        constants.cpp
        timefreq.cpp)

target_link_libraries(chi0_main.exe PRIVATE ${LAPACK_LIBRARIES})

