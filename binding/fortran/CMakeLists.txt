if (BUILD_SHARED_LIBS)
  add_library(rpa_f03_lib SHARED)
else()
  add_library(rpa_f03_lib STATIC)
endif()

set_target_properties(rpa_f03_lib PROPERTIES
  OUTPUT_NAME ${LIBRPA_FORTRAN_LIB_NAME}
  VERSION ${CMAKE_PROJECT_VERSION}
  SOVERSION ${CMAKE_PROJECT_VERSION}
)

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/librpa_f03.f90
  ${CMAKE_CURRENT_BINARY_DIR}/librpa_f03.f90)

target_sources(rpa_f03_lib
  PRIVATE
  ${CMAKE_CURRENT_BINARY_DIR}/librpa_f03.f90
)

target_link_libraries(rpa_f03_lib PRIVATE rpa_lib)

target_include_directories(rpa_f03_lib PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

install(TARGETS rpa_f03_lib
        ARCHIVE
        DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/librpa_f03.mod
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

if(ENABLE_LIBRPA_TEST AND ENABLE_LIBRPA_FORTRAN_TEST)
  add_subdirectory(test)
endif()
