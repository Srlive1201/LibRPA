
if (BUILD_LIBRARY)
  link_libraries(librpa_lib)
endif()

add_executable(test_atoms test_atoms.cpp)
add_executable(test_timefreq test_timefreq.cpp)
add_executable(test_matrix test_matrix.cpp)
add_executable(test_meanfield test_meanfield.cpp)
add_executable(test_complexmatrix test_complexmatrix.cpp)
add_executable(test_pbc test_pbc.cpp)
add_executable(test_parallel_mpi test_parallel_mpi.cpp)
add_executable(test_blas_lapack test_blas_lapack.cpp)

include_directories(PUBLIC ${CMAKE_CURRENT_LIST_DIR})

if (NOT BUILD_LIBRARY)
  target_sources(test_atoms PRIVATE ../src/atoms.cpp)
  target_sources(test_matrix PRIVATE ../src/matrix.cpp)
  target_sources(test_complexmatrix PRIVATE
          ../src/complexmatrix.cpp
          ../src/matrix.cpp
          )
  target_sources(test_meanfield PRIVATE
          ../src/meanfield.cpp
          ../src/constants.cpp
          ../src/complexmatrix.cpp
          ../src/matrix.cpp
          )
  target_sources(test_timefreq PRIVATE
          ../src/constants.cpp
          ../src/timefreq.cpp
          ../src/mathtools.cpp
          ../src/matrix.cpp
          ../src/complexmatrix.cpp
          ${CMAKE_CURRENT_BINARY_DIR}/../src/envs.cpp
          )
  target_sources(test_pbc PRIVATE
          ../src/pbc.cpp
          ../src/vector3_order.cpp
          )
  target_sources(test_blas_lapack PRIVATE
          ../src/complexmatrix.cpp
          ../src/matrix.cpp
          )
  target_sources(test_parallel_mpi PRIVATE
          ../src/parallel_mpi.cpp
          ../src/matrix.cpp
          ../src/complexmatrix.cpp
          ../src/vector3_order.cpp
          )
  target_link_libraries(test_timefreq PRIVATE ${math_libs})
  target_link_libraries(test_matrix PRIVATE ${math_libs})
  target_link_libraries(test_complexmatrix PRIVATE ${math_libs})
  target_link_libraries(test_meanfield PRIVATE ${math_libs})
  target_link_libraries(test_parallel_mpi PRIVATE ${math_libs} ${parallel_libs})
  target_link_libraries(test_blas_lapack PRIVATE ${math_libs})
endif()

add_test(
        NAME test_atoms
        COMMAND $<TARGET_FILE:test_atoms>
)

add_test(
        NAME test_timefreq
        COMMAND $<TARGET_FILE:test_timefreq>
)

add_test(
        NAME test_matrix
        COMMAND $<TARGET_FILE:test_matrix>
)

add_test(
        NAME test_meanfield
        COMMAND $<TARGET_FILE:test_meanfield>
)

add_test(
        NAME test_complexmatrix
        COMMAND $<TARGET_FILE:test_complexmatrix>
)

add_test(
        NAME test_parallel_mpi
        COMMAND mpirun -np 4 $<TARGET_FILE:test_parallel_mpi>
)

add_test(
        NAME test_pbc
        COMMAND $<TARGET_FILE:test_pbc>
)

add_test(
        NAME test_blas_lapack
        COMMAND $<TARGET_FILE:test_blas_lapack>
)

