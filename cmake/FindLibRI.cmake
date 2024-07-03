# Find LibRI headers

# Variables
#   LIBRI_INCLUDE_DIR - specify the path by either environment varaible or -DLIBRI_INCLUDE_DIR to search the cereal headers
#
#   LibRI_FOUND - True if LibRI is found
#   LibRI_INCLUDE_DIR - Where to find the LibRI headers
find_path(LibRI_INCLUDE_DIR
  RI/ri/RI_Tools.h
  HINTS ${LIBRI_INCLUDE_DIR}
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(LibRI
  FOUND_VAR
    LibRI_FOUND
  REQUIRED_VARS
    LibRI_INCLUDE_DIR
  )
