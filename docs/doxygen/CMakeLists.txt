add_doxygen_doc(
  BUILD_DIR
    ${CMAKE_CURRENT_BINARY_DIR}/_build_doxygen
  DOXY_FILE
    ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
  TARGET_NAME
    docs
  COMMENT
    "Doxygen HTML documentation"
  )
