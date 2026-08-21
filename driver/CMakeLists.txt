add_executable(rpa_exe)
set_target_properties(rpa_exe PROPERTIES
  OUTPUT_NAME ${LIBRPA_DRIVER_NAME}
  ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
  LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
  RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})

list(APPEND exe_sources
  input_parser.cpp
  inputfile.cpp
  main.cpp
  task.cpp
  # Specific task implementations
  tasks/task_helper.cpp
  tasks/exx.cpp
  tasks/exx_band.cpp
  tasks/g0w0.cpp
  tasks/g0w0_band.cpp
  tasks/rpa.cpp
  tasks/print_minimax.cpp
  # tasks/test.cpp
  # task_qsgw.cpp
  # task_qsgwA.cpp
  # task_qsgw_band.cpp
  # task_hf_band.cpp
  # task_scRPA.cpp
  # task_scRPA_band.cpp
  # task_screened_coulomb.cpp
  )

target_sources(rpa_exe PRIVATE ${exe_sources})
target_compile_definitions(rpa_exe PRIVATE "LIBRPA_MPI_THREAD_LEVEL=${LIBRPA_MPI_THREAD_LEVEL}")

target_link_libraries(rpa_exe PRIVATE rpa_lib)
if(LIBRPA_BUNDLED_ELPA_PROJECT_TARGET)
  # base_blacs.h includes ELPA headers, so make them available before compiling
  # the driver instead of waiting only for the final rpa_lib link dependency.
  add_dependencies(rpa_exe ${LIBRPA_BUNDLED_ELPA_PROJECT_TARGET})
endif()
if(BUILD_SHARED_LIBS)
  target_link_libraries(rpa_exe PRIVATE ${math_libs} ${parallel_libs})
endif()

install(TARGETS rpa_exe DESTINATION ${CMAKE_INSTALL_BINDIR})

if(LIBRPA_ENABLE_TEST AND LIBRPA_ENABLE_CPP_TEST)
  add_subdirectory(test)
endif()
