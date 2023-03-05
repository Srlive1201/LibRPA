# a simple macro to add source, include, link main library and run test
# see also
# - https://stackoverflow.com/questions/53868015/avoid-repetitive-cmake-code-with-multiple-targets
# - https://stackoverflow.com/questions/19278510/cmake-ctest-add-test-command-appends-extra-quotes-to-my-command-parameters
macro(config_test TARGET NUMPROC)
  add_executable(${TARGET} ${TARGET}.cpp)
  target_link_libraries(${TARGET} PRIVATE rpa_lib)
  if(${NUMPROC} GREATER 1)
    # message("parallel test")
    add_test(
      NAME ${TARGET}
      COMMAND mpirun -np ${NUMPROC} $<TARGET_FILE:${TARGET}>
    )
  else()
    # message("serial test")
    add_test(
      NAME ${TARGET}
      COMMAND $<TARGET_FILE:${TARGET}>
    )
  endif()
endmacro()

# add serial tests
foreach(target
    test_atoms
    test_atomic_basis
    test_timefreq 
    test_matrix
    test_matrix_m
    test_meanfield
    test_complexmatrix
    test_pbc
    test_blas_lapack
    )
  config_test(${target} 1)
endforeach()

# add 4-process tests
foreach(target
  test_parallel_mpi
  test_matrix_m_parallel_utils
  )
  config_test(${target} 4)
endforeach()

# TODO: make use of GTest
