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
       ${CMAKE_CURRENT_SOURCE_DIR}/ekat/cmake/tpls)
  include(EkatBuildEkat)
  set(EKAT_ENABLE_MPI ${HAERO_ENABLE_MPI} CACHE BOOL "Enable MPI")
  set(EKAT_ENABLE_FORTRAN OFF CACHE BOOL "Enable EKAT Fortran support")
  set(EKAT_ENABLE_TESTS OFF CACHE BOOL "Disable EKAT tests")
  set(Kokkos_ENABLE_LIBDL OFF CACHE BOOL "Disable broken Kokkos libdl support")
  set(Kokkos_ENABLE_DEPRECATED_CODE_3 OFF CACHE BOOL "Disable deprecated code to avoid warnings of using deprecated functions.")
  if (HAERO_SKIP_FIND_YAML_CPP)
    set(EKAT_SKIP_FIND_YAML_CPP ON CACHE BOOL "EKAT will build yaml-cpp ")
    message(STATUS "SKIP FIND YAML CPP")
  else()
    set(EKAT_SKIP_FIND_YAML_CPP OFF CACHE BOOL "EKAT will build yaml-cpp ")
    message(STATUS "EKAT WILL LOOK FOR YAML-CPP")
  endif()

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
    if (HAERO_AMD_GPU)
      set(Kokkos_ENABLE_HIP ON CACHE BOOL "Enable HIP Kokkos backend")
    else()
      set(Kokkos_ENABLE_CUDA ON CACHE BOOL "Enable CUDA Kokkos backend")
      set(Kokkos_ENABLE_CUDA_LAMBDA ON CACHE BOOL "Enable CUDA lambdas")
      set(EKAT_NVCC_WRAPPER ${CMAKE_BINARY_DIR}/bin/nvcc_wrapper)
      set(EKAT_NVCC_WRAPPER ${EKAT_NVCC_WRAPPER} PARENT_SCOPE)
      file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/ekat/extern/kokkos/bin/nvcc_wrapper
           DESTINATION ${CMAKE_BINARY_DIR}/bin)
    endif()
    set(Kokkos_ENABLE_SERIAL ON CACHE BOOL "Enable serial Kokkos backend")
  else()
    set(Kokkos_ENABLE_SERIAL ON CACHE BOOL "Enable serial Kokkos backend")
    set(Kokkos_ENABLE_CUDA OFF CACHE BOOL "Disable CUDA Kokkos backend")
    set(Kokkos_ENABLE_HIP OFF CACHE BOOL "Disable HIP Kokkos backend")
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
  if (NOT TARGET yaml-cpp)
    add_library(yaml-cpp STATIC IMPORTED GLOBAL)
    if (CMAKE_BUILD_TYPE MATCHES "Debug")
      set_target_properties(yaml-cpp PROPERTIES
        IMPORTED_LOCATION ${CMAKE_INSTALL_PREFIX}/lib/libyaml-cppd.a)
    else()
      set_target_properties(yaml-cpp PROPERTIES
        IMPORTED_LOCATION ${CMAKE_INSTALL_PREFIX}/lib/libyaml-cpp.a)
    endif()
  endif()
  if (NOT TARGET spdlog)
    add_library(spdlog STATIC IMPORTED GLOBAL)
    if (CMAKE_BUILD_TYPE MATCHES "Debug")
      set_target_properties(spdlog PROPERTIES
        IMPORTED_LOCATION ${CMAKE_INSTALL_PREFIX}/lib/libspdlogd.a)
    else()
      set_target_properties(spdlog PROPERTIES
        IMPORTED_LOCATION ${CMAKE_INSTALL_PREFIX}/lib/libspdlog.a)
    endif()
  endif()
endif()

add_library(ext_libraries INTERFACE)
target_link_libraries (ext_libraries INTERFACE ekat;${HAERO_EXT_LIBRARIES})

# Add all the libraries to the list maintained for the standalone build config file.
set(HAERO_LIBRARIES ekat;${HAERO_EXT_LIBRARIES};${HAERO_LIBRARIES} PARENT_SCOPE)
list(REMOVE_DUPLICATES HAERO_EXT_INCLUDE_DIRS)
set(HAERO_EXT_INCLUDE_DIRS ${HAERO_EXT_INCLUDE_DIRS} PARENT_SCOPE)

#----------------------
# Installation targets
#----------------------

install(DIRECTORY ${PROJECT_BINARY_DIR}/include/ DESTINATION include)
install(DIRECTORY ${PROJECT_BINARY_DIR}/bin/ DESTINATION bin)
