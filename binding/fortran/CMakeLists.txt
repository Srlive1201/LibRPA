if (BUILD_LIBRPA_SHARED)
  add_library(rpa_f_lib SHARED)
else()
  add_library(rpa_f_lib STATIC)
endif()

set_target_properties(rpa_f_lib PROPERTIES
  OUTPUT_NAME ${FORTRAN_LIB_NAME}
  VERSION ${CMAKE_PROJECT_VERSION}
  SOVERSION ${CMAKE_PROJECT_VERSION}
)

target_sources(rpa_f_lib
  PRIVATE
  librpa_f.f90
)
target_link_libraries(rpa_f_lib PRIVATE rpa_lib)

install(TARGETS rpa_f_lib
        ARCHIVE
        DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/librpa_f.mod
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
