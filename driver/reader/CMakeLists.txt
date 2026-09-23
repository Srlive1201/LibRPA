add_library(librpa_file_reader
  file_reader.cpp
  read_data.cpp
  reader_band.cpp
  reader_structure.cpp
  reader_basis.cpp
  reader_lri.cpp
  reader_coulomb.cpp
  reader_eigenvec.cpp)
target_link_libraries(librpa_file_reader PUBLIC rpa_lib)
if(LIBRPA_BUNDLED_ELPA_PROJECT_TARGET)
  # Reader headers include base_blacs.h, which needs the generated ELPA headers.
  add_dependencies(librpa_file_reader ${LIBRPA_BUNDLED_ELPA_PROJECT_TARGET})
endif()
set_target_properties(librpa_file_reader PROPERTIES OUTPUT_NAME rpa_file_reader)
install(TARGETS librpa_file_reader
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
