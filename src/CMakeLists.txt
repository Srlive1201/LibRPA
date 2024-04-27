# By default, if the compilation succeeds, you will have the executable at the "src" subdirectory in the build directory
# To put it in the build root directory
# set_target_properties(rpa_exe
#   PROPERTIES
#   ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
#   LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
#   RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})

if (BUILD_LIBRPA_SHARED)
  add_library(rpa_lib SHARED)
else()
  add_library(rpa_lib STATIC)
endif()

set_target_properties(rpa_lib PROPERTIES
  OUTPUT_NAME ${LIB_NAME}
  VERSION ${CMAKE_PROJECT_VERSION}
  SOVERSION ${CMAKE_PROJECT_VERSION}
)

configure_file(envs.cpp.in envs.cpp)
list(APPEND lib_sources
        # aperiodic_chi0.cpp
        analycont.cpp
        atomic_basis.cpp
        # cal_periodic_chi0.cpp
        chi0.cpp
        complexmatrix.cpp
        constants.cpp
        coulmat.cpp
        dielecmodel.cpp
        ${CMAKE_CURRENT_BINARY_DIR}/envs.cpp
        epsilon.cpp
        exx.cpp
        fitting.cpp
        gw.cpp
        librpa.cpp
        librpa_main.cpp
        interpolate.cpp
        mathtools.cpp
        matrix.cpp
        matrix3.cpp
        meanfield.cpp
        parallel_mpi.cpp
        params.cpp
        pbc.cpp
        profiler.cpp
        qpe_solver.cpp
        ri.cpp
        timefreq.cpp
        vector3_order.cpp
        )

if(USE_GREENX_MINIMAX)
  list(APPEND lib_sources get_minimax_gx.cpp gx_minimax_wrp.f90)
  target_link_libraries(rpa_lib PRIVATE GXCommon LibGXMiniMax)
  target_include_directories(rpa_lib
    PRIVATE ${PROJECT_BINARY_DIR}/thirdparty/greenX-6ff8a00/modules
    )
else()
  list(APPEND lib_sources get_minimax_local.cpp)
endif()

target_include_directories(rpa_lib
  PRIVATE
    # ${CMAKE_CURRENT_LIST_DIR}
    ${PROJECT_SOURCE_DIR}/src
  )

if(USE_LIBRI)
  target_include_directories(rpa_lib
    PRIVATE
      ${Cereal_INCLUDE_DIR}
      ${LibComm_INCLUDE_DIR}
      ${LibRI_INCLUDE_DIR}
    )
endif()

# source must be PRIVATE to avoid recompilation when linking driver to this library
target_sources(rpa_lib PRIVATE ${lib_sources})
target_link_libraries(rpa_lib PUBLIC ${math_libs} ${parallel_libs})

install(TARGETS rpa_lib
        ARCHIVE
        DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(FILES ${PROJECT_SOURCE_DIR}/src/librpa.h
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

if(ENABLE_TEST)
  add_subdirectory(test)
endif()
