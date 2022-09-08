# Find Cereal headers
# Adapted from the ABACUS FindCereal module and cmake cookbook ch03-r10
#
# Variables
#   CEREAL_INCLUDE_DIR - specify the path by either environment varaible or -DCEREAL_INCLUDE_DIR to search the cereal headers
#
#   Cereal_FOUND - True if cereal is found.
#   Cereal_INCLUDE_DIR - Where to find cereal headers.

# check environment variable if cmake definition is not set
if(NOT CEREAL_INCLUDE_DIR)
  set(CEREAL_INCLUDE_DIR $ENV{CEREAL_INCLUDE_DIR})
endif()

find_path(Cereal_INCLUDE_DIR
  cereal/cereal.hpp
  HINTS ${CEREAL_INCLUDE_DIR}
)

# Disable auto download. Handle CEREAL_FOUND instead
# if(NOT Cereal_INCLUDE_DIR)
#     include(FetchContent)
#     FetchContent_Declare(
#         cereal
#         GIT_REPOSITORY https://github.com.cnpmjs.org/USCiLab/cereal.git
#         # URL https://github.com/USCiLab/cereal/archive/refs/tags/v1.3.0.tar.gz
#     )
#     set(Cereal_INCLUDE_DIR ${cereal_SOURCE_DIR})
# endif()

# Handle the QUIET and REQUIRED arguments and
# set Cereal_FOUND to TRUE if all variables are non-zero.
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Cereal
  FOUND_VAR
    Cereal_FOUND
  REQUIRED_VARS
    Cereal_INCLUDE_DIR
  )

# Copy the results to the output variables and target.
# mark_as_advanced(Cereal_INCLUDE_DIR)
