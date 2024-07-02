add_executable(rpa_exe)
set_target_properties(rpa_exe PROPERTIES
  OUTPUT_NAME ${DRIVER_NAME}
  ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
  LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
  RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})

list(APPEND exe_sources
  main.cpp
  read_data.cpp
  write_aims.cpp
  inputfile.cpp
  task_rpa.cpp
  )

target_include_directories(rpa_exe
  PUBLIC
    ${PROJECT_SOURCE_DIR}/src
  )

target_sources(rpa_exe PRIVATE ${exe_sources})

target_link_libraries(rpa_exe PRIVATE rpa_lib ${math_libs} ${parallel_libs})

install(TARGETS rpa_exe DESTINATION ${CMAKE_INSTALL_BINDIR})
