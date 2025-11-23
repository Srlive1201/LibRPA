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

configure_file(envs_cmake.cpp.in envs_cmake.cpp)
list(APPEND lib_sources
        # aperiodic_chi0.cpp
        analycont.cpp

        api/app_exx.cpp
        api/app_rpa.cpp
        api/app_gw.cpp
        api/librpa_compute.cpp
        api/librpa_global.cpp
        api/librpa_handler.cpp
        api/librpa_options.cpp
        api/librpa.cpp
        api/librpa_main.cpp

        atomic_basis.cpp
        # cal_periodic_chi0.cpp
        chi0.cpp
        constants.cpp
        coulmat.cpp
        dielecmodel.cpp
        ${CMAKE_CURRENT_BINARY_DIR}/envs_cmake.cpp

        global/instance_manager.cpp

        math/complexmatrix.cpp
        math/fitting.cpp
        math/interpolate.cpp
        math/mathtools.cpp
        math/matrix.cpp
        math/matrix3.cpp
        math/utils_matrix_mpi.cpp
        math/vector3_order.cpp

        mpi/base_blacs.cpp
        mpi/base_mpi.cpp
        mpi/envs_blacs.cpp
        mpi/envs_mpi.cpp
        mpi/utils_blacs.cpp

        envs_io.cpp
        epsilon.cpp
        exx.cpp
        geometry.cpp
        gw.cpp
        meanfield.cpp
        params.cpp
        pbc.cpp
        profiler.cpp
        qpe_solver.cpp
        ri.cpp
        timefreq.cpp
        utils_atomic_basis_blacs.cpp
        utils_cmake.cpp
        utils_mem.cpp
        utils_timefreq.cpp
)

if(USE_GREENX_API)
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
    ${PROJECT_SOURCE_DIR}/include
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
install(FILES
        ${PROJECT_SOURCE_DIR}/src/api/librpa.h
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

if(ENABLE_TEST)
  add_subdirectory(test)
endif()
