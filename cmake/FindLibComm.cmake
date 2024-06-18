# Find LibComm headers

# Variables
#   LIBCOMM_INCLUDE_DIR - specify the path by -DLIBCOMM_INCLUDE_DIR to search the LibComm headers
#
# Returns
#   LibComm_FOUND - True if LibComm is found
#   LibComm_INCLUDE_DIR - Where to find the LibComm headers

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
