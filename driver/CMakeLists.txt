add_executable(${PROJECT_NAME}_exe)
set_target_properties(${PROJECT_NAME}_exe PROPERTIES
  OUTPUT_NAME ${LIBRPA_DRIVER_NAME}
  ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
  LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
  RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})

list(APPEND exe_sources
  driver.cpp
  inputfile.cpp
  main.cpp
  read_data.cpp
  task.cpp
  # Specific task implementations
  task/exx.cpp
  task/exx_band.cpp
  task/g0w0.cpp
  task/gw_band.cpp
  task/rpa.cpp
  # task/screened_coulomb.cpp
  task/print_minimax.cpp
  # task/test.cpp
)

target_include_directories(${PROJECT_NAME}_exe PRIVATE ${PROJECT_NAME}_lib)
target_sources(${PROJECT_NAME}_exe PRIVATE ${exe_sources})

target_link_libraries(${PROJECT_NAME}_exe PRIVATE ${PROJECT_NAME}_lib)
if(BUILD_SHARED_LIBS)
  target_link_libraries(${PROJECT_NAME}_exe PRIVATE ${math_libs} ${parallel_libs})
endif()

install(TARGETS ${PROJECT_NAME}_exe DESTINATION ${CMAKE_INSTALL_BINDIR})
