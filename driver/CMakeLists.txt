add_executable(rpa_exe)
set_target_properties(rpa_exe PROPERTIES
  OUTPUT_NAME ${DRIVER_NAME}
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
  # task/exx.cpp
  # task/exx_band.cpp
  task/gw.cpp
  # task/gw_band.cpp
  # task/rpa.cpp
  # task/screened_coulomb.cpp
  task/print_minimax.cpp
  # task/test.cpp
  )

target_include_directories(rpa_exe
  PUBLIC
    ${PROJECT_SOURCE_DIR}/src
  )

target_sources(rpa_exe PRIVATE ${exe_sources})
if(USE_LIBRI)
  target_include_directories(rpa_exe
    PRIVATE
      ${Cereal_INCLUDE_DIR}
      ${LibComm_INCLUDE_DIR}
      ${LibRI_INCLUDE_DIR}
    )
endif()

target_link_libraries(rpa_exe PRIVATE rpa_lib ${math_libs} ${parallel_libs})

install(TARGETS rpa_exe DESTINATION ${CMAKE_INSTALL_BINDIR})
