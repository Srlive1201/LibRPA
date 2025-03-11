if (BUILD_LIBRPA_SHARED)
  add_library(qsgw_lib SHARED)
else()
  add_library(qsgw_lib STATIC)
endif()

set_target_properties(qsgw_lib PROPERTIES
  OUTPUT_NAME qsgw
  VERSION ${CMAKE_PROJECT_VERSION}
  SOVERSION ${CMAKE_PROJECT_VERSION}
)

list(APPEND qsgw_sources
  Hamiltonian.cpp
  convert_csc.cpp
  fermi_energy_occupation.cpp
)

target_include_directories(qsgw_lib
  PRIVATE
    ${PROJECT_SOURCE_DIR}/src
    ${PROJECT_SOURCE_DIR}/qsgw
  )

target_sources(qsgw_lib PRIVATE ${qsgw_sources})
target_link_libraries(qsgw_lib PRIVATE rpa_lib ${math_libs} ${parallel_libs})
