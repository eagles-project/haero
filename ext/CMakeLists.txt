include(ExternalProject)
include(GNUInstallDirs)
include_directories(${PROJECT_BINARY_DIR}/include)

# Create canonical top-level directories in the build directory.
foreach(dir include;lib;bin;share)
  set(top_level_dir ${PROJECT_BINARY_DIR}/${dir})
  if (NOT EXISTS ${top_level_dir})
    execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${top_level_dir})
  endif()
endforeach()

#----------------------------------------------------------------------------
# E3SM Kokkos Application Toolkit (EKAT) library and friends.
#----------------------------------------------------------------------------
if (HAERO_BUILDS_EKAT)
  list(APPEND CMAKE_MODULE_PATH
       ${CMAKE_CURRENT_SOURCE_DIR}/ekat/cmake
       ${CMAKE_CURRENT_SOURCE_DIR}/ekat/cmake/pkg_build)
  include(EkatBuildEkat)
  set(EKAT_ENABLE_MPI ${HAERO_ENABLE_MPI} CACHE BOOL "Enable MPI")
  set(EKAT_ENABLE_TESTS OFF CACHE BOOL "Disable EKAT tests")
  set(Kokkos_ENABLE_LIBDL OFF CACHE BOOL "Disable broken Kokkos libdl support")
  set(Kokkos_ENABLE_DEPRECATED_CODE_3 OFF CACHE BOOL "Disable deprecated code to avoid warnings of using deprecated functions.")
  if (CMAKE_BUILD_TYPE STREQUAL Release)
    set(Kokkos_ENABLE_DEBUG FALSE)
    set(Kokkos_ENABLE_AGGRESSIVE_VECTORIZATION ON)
  else()
    set(Kokkos_ENABLE_DEBUG TRUE)
    set(Kokkos_ENABLE_AGGRESSIVE_VECTORIZATION OFF)
  endif()
  if (APPLE)
    # floating point exceptions don't work properly on Macs
    set(EKAT_ENABLE_FPE OFF CACHE BOOL "")
  endif()
  if (HAERO_ENABLE_GPU)
    set(Kokkos_ENABLE_SERIAL ON CACHE BOOL "Enable serial Kokkos backend")
    set(Kokkos_ENABLE_CUDA ON CACHE BOOL "Enable CUDA Kokkos backend")
    set(Kokkos_ENABLE_CUDA_LAMBDA ON CACHE BOOL "Enable CUDA lambdas")
  else()
    set(Kokkos_ENABLE_SERIAL ON CACHE BOOL "Enable serial Kokkos backend")
    set(Kokkos_ENABLE_CUDA OFF CACHE BOOL "Disable CUDA Kokkos backend")
    if (OPENMP_FOUND)
      message(STATUS "OpenMP is enabled")
      set(Kokkos_ENABLE_OPENMP ON CACHE BOOL "Enable OpenMP Kokkos backend")
    else()
      message(STATUS "OpenMP is not supported. Using threads backend instead.")
      set(Kokkos_ENABLE_OPENMP OFF CACHE BOOL "Disable OpenMP Kokkos backend")
      set(Kokkos_ENABLE_PTHREAD ON CACHE BOOL "Enable pthreads Kokkos backend")
    endif()
  endif()
  BuildEkat(PREFIX "HAERO")
  list(APPEND EKAT_INCLUDE_DIRS
       ${CMAKE_CURRENT_SOURCE_DIR}/ekat/src
       ${PROJECT_BINARY_DIR}/externals/ekat/src
       ${CMAKE_CURRENT_SOURCE_DIR}/ekat/extern/kokkos/tpls/desul/include
       ${CMAKE_CURRENT_SOURCE_DIR}/ekat/extern/kokkos/core/src
       ${CMAKE_CURRENT_SOURCE_DIR}/ekat/extern/kokkos/containers/src
       ${PROJECT_BINARY_DIR}/externals/kokkos
       ${CMAKE_CURRENT_SOURCE_DIR}/ekat/extern/spdlog/include
       ${CMAKE_CURRENT_SOURCE_DIR}/ekat/extern/yaml-cpp/include)
  file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/ekat/extern/kokkos/bin/nvcc_wrapper
       DESTINATION ${CMAKE_BINARY_DIR}/bin)
else()
  # FIXME: this is horrible. We need to package EKAT better so it includes
  # FIXME: the various information and files that we cobble together here.

  # Copy EKAT's test-launcher and nvcc_wrapper into place.
  file(COPY ${EKAT_SOURCE_DIR}/bin/test-launcher DESTINATION ${CMAKE_BINARY_DIR}/bin)
  if (HAERO_ENABLE_GPU)
    file(COPY ${EKAT_SOURCE_DIR}/bin/nvcc_wrapper DESTINATION ${CMAKE_BINARY_DIR}/bin)
  endif()
  # Copy Catch2 single header to where Haero can reach it.
  file(COPY ${EKAT_SOURCE_DIR}/extern/Catch2/single_include/catch2 DESTINATION ${PROJECT_BINARY_DIR}/include)
  # Copy Kokkos_Random.hpp header to where Haero can reach it.
  file(COPY
       ${EKAT_SOURCE_DIR}/extern/kokkos/algorithms/src/Kokkos_Random.hpp
       ${EKAT_SOURCE_DIR}/extern/kokkos/algorithms/src/Kokkos_Sort.hpp
       DESTINATION ${PROJECT_BINARY_DIR}/include)
  # Copy EKAT cmake functions to where Haero can see them.
  foreach(script EkatCreateUnitTest EkatUtils)
    file(COPY ${EKAT_SOURCE_DIR}/cmake/${script}.cmake DESTINATION ${PROJECT_BINARY_DIR}/share)
  endforeach()
  list(APPEND EKAT_INCLUDE_DIRS
       ${EKAT_SOURCE_DIR}/src
       ${EKAT_BINARY_DIR}/src
       ${EKAT_SOURCE_DIR}/extern/kokkos/tpls/desul/include
       ${EKAT_SOURCE_DIR}/extern/kokkos/core/src
       ${EKAT_SOURCE_DIR}/extern/kokkos/containers/src
       ${EKAT_BINARY_DIR}/../kokkos
       ${EKAT_SOURCE_DIR}/extern/spdlog/include
       ${EKAT_SOURCE_DIR}/extern/yaml-cpp/include
       ${EKAT_BINARY_DIR}/extern/Catch2/single_include)
  set(EKAT_TPL_LIBRARIES kokkos;yaml-cpp;spdlog)
  add_library(ekat STATIC IMPORTED GLOBAL)
  set_target_properties(ekat PROPERTIES
    IMPORTED_LOCATION ${EKAT_BINARY_DIR}/src/ekat/libekat.a)
  add_library(ekat_test_main STATIC IMPORTED GLOBAL)
  set_target_properties(ekat_test_main PROPERTIES
    IMPORTED_LOCATION ${EKAT_BINARY_DIR}/src/ekat/libekat_test_main.a)
  add_library(ekat_test_session STATIC IMPORTED GLOBAL)
  set_target_properties(ekat_test_session PROPERTIES
    IMPORTED_LOCATION ${EKAT_BINARY_DIR}/src/ekat/libekat_test_session.a)

  add_library(yaml-cpp STATIC IMPORTED GLOBAL)
  add_library(spdlog STATIC IMPORTED GLOBAL)
  if (CMAKE_BUILD_TYPE MATCHES "Debug")
    set_target_properties(yaml-cpp PROPERTIES
      IMPORTED_LOCATION ${EKAT_BINARY_DIR}/../yaml-cpp/libyaml-cppd.a)
    set_target_properties(spdlog PROPERTIES
      IMPORTED_LOCATION ${EKAT_BINARY_DIR}/../spdlog/libspdlogd.a)
  else()
    set_target_properties(yaml-cpp PROPERTIES
      IMPORTED_LOCATION ${EKAT_BINARY_DIR}/../yaml-cpp/libyaml-cpp.a)
    set_target_properties(spdlog PROPERTIES
      IMPORTED_LOCATION ${EKAT_BINARY_DIR}/../spdlog/libspdlog.a)
  endif()
  add_library(kokkoscore STATIC IMPORTED GLOBAL)
  set_target_properties(kokkoscore PROPERTIES
      IMPORTED_LOCATION ${EKAT_BINARY_DIR}/../kokkos/core/src/libkokkoscore.a)
  add_library(kokkoscontainers STATIC IMPORTED GLOBAL)
  set_target_properties(kokkoscontainers PROPERTIES
      IMPORTED_LOCATION ${EKAT_BINARY_DIR}/../kokkos/containers/src/libkokkoscontainers.a)
  add_library(kokkossimd STATIC IMPORTED GLOBAL)
  set_target_properties(kokkossimd PROPERTIES
      IMPORTED_LOCATION ${EKAT_BINARY_DIR}/../kokkos/simd/src/libkokkossimd.a)

  # set up the "kokkos" interface library
  add_library(kokkos INTERFACE)
  add_library(Kokkos::kokkos ALIAS kokkos)
  target_link_libraries(kokkos INTERFACE kokkoscore kokkoscontainers)
endif()
set(EKAT_NVCC_WRAPPER ${CMAKE_BINARY_DIR}/bin/nvcc_wrapper)
set(EKAT_NVCC_WRAPPER ${EKAT_NVCC_WRAPPER} PARENT_SCOPE)

set(HAERO_EXT_LIBRARIES ekat;${EKAT_TPL_LIBRARIES};${HAERO_EXT_LIBRARIES})
foreach(inc_dir ${EKAT_INCLUDE_DIRS})
  list(APPEND HAERO_EXT_INCLUDE_DIRS ${inc_dir})
endforeach()

# This target can be used as a dependency to ensure all above libraries are
# built before proceeding.
add_custom_target(ext_libraries DEPENDS ${HAERO_EXT_LIBRARIES})

# Add all the libraries to the build system at large.
set(HAERO_LIBRARIES ${HAERO_EXT_LIBRARIES};${HAERO_LIBRARIES} PARENT_SCOPE)
list(REMOVE_DUPLICATES HAERO_EXT_INCLUDE_DIRS)
set(HAERO_EXT_INCLUDE_DIRS ${HAERO_EXT_INCLUDE_DIRS} PARENT_SCOPE)

#----------------------
# Installation targets
#----------------------

install(DIRECTORY ${PROJECT_BINARY_DIR}/include/ DESTINATION include)
install(DIRECTORY ${PROJECT_BINARY_DIR}/bin/ DESTINATION bin)

