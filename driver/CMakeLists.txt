add_executable(rpa_exe)
set_target_properties(rpa_exe PROPERTIES
  OUTPUT_NAME ${LIBRPA_DRIVER_NAME}
  ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
  LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
  RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})

list(APPEND exe_sources
  driver.cpp
  inputfile.cpp
  main.cpp
  read_data.cpp
  write_aims.cpp
  task.cpp
  # Specific task implementations
  task/exx.cpp
  # task/exx_band.cpp
  task/g0w0.cpp
  # task/gw_band.cpp
  task/rpa.cpp
  # task/screened_coulomb.cpp
  task/print_minimax.cpp
  # task/test.cpp
  )

target_include_directories(rpa_exe PRIVATE rpa_lib)
target_sources(rpa_exe PRIVATE ${exe_sources})

target_link_libraries(rpa_exe PRIVATE rpa_lib)
if(BUILD_SHARED_LIBS)
  target_link_libraries(rpa_exe PRIVATE ${math_libs} ${parallel_libs})
endif()

install(TARGETS rpa_exe DESTINATION ${CMAKE_INSTALL_BINDIR})
