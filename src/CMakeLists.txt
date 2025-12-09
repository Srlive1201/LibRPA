# By default, if the compilation succeeds, you will have the executable at the "src" subdirectory in the build directory
# To put it in the build root directory
# set_target_properties(rpa_exe
#   PROPERTIES
#   ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
#   LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
#   RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})

if (BUILD_SHARED_LIBS)
  add_library(rpa_lib SHARED)
else()
  add_library(rpa_lib STATIC)
endif()

add_dependencies(rpa_lib LibGXMiniMax)

set_target_properties(rpa_lib PROPERTIES
  OUTPUT_NAME ${LIBRPA_LIB_NAME}
  VERSION ${CMAKE_PROJECT_VERSION}
  SOVERSION ${CMAKE_PROJECT_VERSION}
)

target_link_libraries(rpa_lib PRIVATE LibGXMiniMax)
target_include_directories(rpa_lib PRIVATE
  $<TARGET_PROPERTY:LibGXMiniMax,Fortran_MODULE_DIRECTORY>
)

target_include_directories(rpa_lib
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${PROJECT_SOURCE_DIR}/include
  )
configure_file(utils/envs_cmake.cpp.in envs_cmake.cpp)
configure_file(utils/envs_cmake.h envs_cmake.h)  # only to make the configured envs_cmake.cpp compile

list(APPEND lib_sources
        ${CMAKE_CURRENT_BINARY_DIR}/envs_cmake.cpp

        # api/app_exx.cpp
        # api/app_rpa.cpp
        # api/app_gw.cpp
        api/compute_exx.cpp
        api/compute_g0w0.cpp
        api/compute_rpa.cpp
        api/dataset.cpp
        api/dataset_helper.cpp
        api/global.cpp
        api/handler.cpp
        api/input.cpp
        api/instance_manager.cpp
        api/librpa.cpp
        api/options.cpp

        core/analycont.cpp
        core/atomic_basis.cpp
        core/chi0.cpp
        core/coulmat.cpp
        core/dielecmodel.cpp
        core/epsilon.cpp
        core/exx.cpp
        core/geometry.cpp
        core/gw.cpp
        core/meanfield.cpp
        core/meanfield_mpi.cpp
        core/minimax_gx.cpp
        core/pbc.cpp
        core/qpe_solver.cpp
        core/ri.cpp
        core/timefreq.cpp
        core/utils_atomic_basis_blacs.cpp

        io/global_io.cpp
        io/fs.cpp

        interface/gx_minimax_wrp.f90

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
        mpi/utils_blacs.cpp
        mpi/global_mpi.cpp

        utils/constants.cpp
        utils/profiler.cpp
        utils/utils_cmake.cpp
        utils/utils_mem.cpp
)

if(LIBRPA_USE_LIBRI)
  target_include_directories(rpa_lib
    PUBLIC
      ${Cereal_INCLUDE_DIR}
      ${LibComm_INCLUDE_DIR}
      ${LibRI_INCLUDE_DIR}
    )
endif()

# source must be PRIVATE to avoid recompilation when linking driver to this library
target_sources(rpa_lib PRIVATE ${lib_sources})
target_link_libraries(rpa_lib PRIVATE ${math_libs} ${parallel_libs})

install(TARGETS rpa_lib
        ARCHIVE
        DESTINATION ${CMAKE_INSTALL_LIBDIR})

if(ENABLE_LIBRPA_TEST AND ENABLE_LIBRPA_CPP_TEST)
  add_subdirectory(test)
endif()
