add_executable(test_atoms test_atoms.cpp)
add_executable(test_timefreq test_timefreq.cpp)
add_executable(test_matrix test_matrix.cpp)

include_directories(PUBLIC ${CMAKE_CURRENT_LIST_DIR})
# target_include_directories(test_atoms
#   PUBLIC
#     ${CMAKE_CURRENT_LIST_DIR})

if (BUILD_LIBRARY)
  # target_link_libraries(test_atoms librpa_lib)
  link_libraries(librpa_lib)
else()
  # target_link_libraries(test_atoms ${LAPACK_LIBRARIES})
  target_sources(test_atoms PRIVATE ../src/atoms.cpp)
  target_sources(test_matrix PRIVATE ../src/matrix.cpp)
  target_sources(test_timefreq PRIVATE 
    ../src/timefreq.cpp
    ../src/matrix.cpp
    ../src/complexmatrix.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/../src/envs.cpp
    )
  target_link_libraries(test_timefreq PRIVATE ${LAPACK_LIBRARIES})
  target_link_libraries(test_matrix PRIVATE ${LAPACK_LIBRARIES})
endif()

add_test(
  NAME test_atoms
  COMMAND $<TARGET_FILE:test_atoms>
  )

add_test(
  NAME test_timefreq
  COMMAND $<TARGET_FILE:test_timefreq>
  )

add_test(
  NAME test_matrix
  COMMAND $<TARGET_FILE:test_matrix>
  )

