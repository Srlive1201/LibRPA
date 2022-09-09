# Find LibComm headers

# Variables
#   LIBCOMM_INCLUDE_DIR - specify the path by either environment varaible or -DLIBCOMM_INCLUDE_DIR to search the cereal headers
#
#   LibComm_FOUND - True if LibComm is found
#   LibComm_INCLUDE_DIR - Where to find the LibComm headers

# check environment variable if cmake definition is not set
if(NOT LIBCOMM_INCLUDE_DIR)
  set(LIBCOMM_INCLUDE_DIR $ENV{LIBCOMM_INCLUDE_DIR})
endif()

find_path(LibComm_INCLUDE_DIR
  Comm/Comm_Tools.h
  HINTS ${LIBCOMM_INCLUDE_DIR}
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(LibComm
  FOUND_VAR
    LibComm_FOUND
  REQUIRED_VARS
    LibComm_INCLUDE_DIR
  )
