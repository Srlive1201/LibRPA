add_executable(
  test_input_parser
  test_input_parser.cpp
  ../input_parser.cpp
)
add_test(NAME test_input_parser COMMAND $<TARGET_FILE:test_input_parser>)

add_executable(test_driver_params test_driver_params.cpp)
target_link_libraries(test_driver_params PRIVATE rpa_lib)
add_test(NAME test_driver_params COMMAND $<TARGET_FILE:test_driver_params>)

add_executable(test_velocity_reader test_velocity_reader.cpp)
target_link_libraries(test_velocity_reader PRIVATE rpa_lib)
add_test(NAME test_velocity_reader COMMAND $<TARGET_FILE:test_velocity_reader>)
