macro(config_test_f03 TARGET NUMPROC)
  add_executable(${TARGET} ${TARGET}.f90)
  add_dependencies(${TARGET} rpa_f03_lib rpa_lib)
  target_link_libraries(${TARGET} PRIVATE rpa_f03_lib rpa_lib)
  if(BUILD_SHARED_LIBS)
    target_link_libraries(${TARGET} PRIVATE ${math_libs} ${parallel_libs})
  endif()
  target_include_directories(test_f03_binding PRIVATE
    $<TARGET_PROPERTY:rpa_f03_lib,Fortran_MODULE_DIRECTORY>
  )
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


config_test_f03(test_f03_binding 1)
