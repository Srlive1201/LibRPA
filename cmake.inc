# Basic C++ compiler setup
# set(CMAKE_C_COMPILER "mpiicc" CACHE STRING "" FORCE)
set(CMAKE_CXX_COMPILER "mpiicpx" CACHE STRING "" FORCE)
set(CMAKE_Fortran_COMPILER "ifx" CACHE STRING "" FORCE)




set(CMAKE_CXX_FLAGS "-g -O2 -fopenmp -DLIBRPA_VERBOSE" CACHE STRING "" FORCE)
if(DEBUG)
  set(CMAKE_CXX_FLAGS "-g -O0 -Wall -fopenmp" CACHE STRING "" FORCE)
endif()

set(CEREAL_INCLUDE_DIR "/home/linpz/software/cereal/cereal-1.3.0/include")
set(LIBRI_INCLUDE_DIR "/home/ghj/abacus/GW_240302/LibRI-loop3/include")
# set(LIBRI_INCLUDE_DIR "/home/ghj/abacus/LibRI_time-out/include")
set(LIBCOMM_INCLUDE_DIR "/home/ghj/abacus/241111_librpa/LibComm/include")



set(CEREAL_INCLUDE_DIR "/home/elhacedor/LibRPA/thirdparty/cereal-1.3.0/include")
# set(LIBRI_INCLUDE_DIR "/home/zhangmy/programs/LibRI-loop3-mygh/include")
# set(LIBRI_INCLUDE_DIR "/home/bhjia/LibRI-wait/include")

set(LIBRI_INCLUDE_DIR "/home/elhacedor/LibRPA/thirdparty/LibRI/include")

set(LIBCOMM_INCLUDE_DIR "/home/elhacedor/LibRPA/thirdparty/LibComm/include")
# set(LAPACK_LIBRARIES "-L$ENV{MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core")
set(LAPACK_LIBRARIES "-L$ENV{MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lmkl_blacs_intelmpi_lp64 -lpthread -lm -ldl")

# whether to use LibRI for tensor contraction
# set(USE_LIBRI ON)
