# extract the executable name and library name.
file(READ "${PROJECT_SOURCE_DIR}/src/Makefile" mkfile)
string(REGEX MATCH " *EXENAME *= *[0-9a-zA-Z._-]*" exename ${mkfile})
string(REGEX REPLACE " *EXENAME *= *" "" exename ${exename})
string(REGEX MATCH " *LIBNAME *= *[0-9a-zA-Z._-]*" libname ${mkfile})
string(REGEX REPLACE " *LIBNAME *= *" "" libname ${libname})

# add the executable target, with the main function code
add_executable(librpa_exe lib_main.cpp)
set_target_properties(librpa_exe PROPERTIES
  OUTPUT_NAME ${exename})
# Specify where the executable will be generated.
# By default, if the compilation succeeds, you will have the executable at the "src" subdirectory in the build directory
# To put it in the build root directory
set_target_properties(librpa_exe
  PROPERTIES
  ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
  LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
  RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})

target_include_directories(librpa_exe
  PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}
  )

if(USE_LIBRI)
  target_include_directories(librpa_exe
    PUBLIC
      ${Cereal_INCLUDE_DIR}
      ${LibComm_INCLUDE_DIR}
      ${LibRI_INCLUDE_DIR}
    )
endif()

configure_file(envs.cpp.in envs.cpp)
list(APPEND sources
        aperiodic_chi0.cpp
        cal_periodic_chi0.cpp
        chi0.cpp
        complexmatrix.cpp
        constants.cpp
        coulmat.cpp
        ${CMAKE_CURRENT_BINARY_DIR}/envs.cpp
        epsilon.cpp
        exx.cpp
        Gauss_Quadrature.cpp
        input.cpp
        mathtools.cpp
        matrix.cpp
        matrix3.cpp
        meanfield.cpp
        parallel_mpi.cpp
        params.cpp
        pbc.cpp
        profiler.cpp
        read_aims.cpp
        ri.cpp
        timefreq.cpp
        vector3_order.cpp
	d_chi0.cpp
        )


if(BUILD_LIBRARY)
  add_library(librpa_lib)
  set_target_properties(librpa_lib PROPERTIES
    OUTPUT_NAME ${libname}
    ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
    LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
    RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
  target_include_directories(librpa_lib
    PUBLIC
      ${CMAKE_CURRENT_LIST_DIR}
    )
  if(USE_LIBRI)
    target_include_directories(librpa_lib
      PUBLIC
        ${Cereal_INCLUDE_DIR}
        ${LibComm_INCLUDE_DIR}
        ${LibRI_INCLUDE_DIR}
      )
  endif()
  target_sources(librpa_lib PUBLIC ${sources})
  target_link_libraries(librpa_lib PUBLIC ${math_libs} ${parallel_libs})
  target_link_libraries(librpa_exe PRIVATE librpa_lib)
else()
  target_sources(librpa_exe PRIVATE ${sources})
  target_link_libraries(librpa_exe PRIVATE ${math_libs} ${parallel_libs})
endif()

