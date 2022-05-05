# add_executable(test_pbc test_pbc.cpp)
#
# target_include_directories(test_pbc
#   PUBLIC
#     ${CMAKE_CURRENT_LIST_DIR})
# if (BUILD_LIBRARY)
#   target_link_libraries(test_pbc cmake_cpp_example_lib)
# else()
#   # target_sources(test_pbc PRIVATE ${PROJECT_SOURCE_DIR}/src/pbc.cpp)
#   target_link_libraries(test_pbc ${LAPACK_LIBRARIES})
# endif()
#
# add_test(
#   NAME test_pbc
#   COMMAND $<TARGET_FILE:test_pbc>
#   )
