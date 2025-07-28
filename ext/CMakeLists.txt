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

  # Set Kokkos config options
  set(Kokkos_ENABLE_LIBDL OFF CACHE BOOL "Whether to enable the LIBDL library")
  set(Kokkos_ENABLE_DEPRECATED_CODE_4 OFF CACHE BOOL "Whether code deprecated in major release 4 is available")

  if (CMAKE_BUILD_TYPE STREQUAL Release)
    set(Kokkos_ENABLE_DEBUG FALSE)
    set(Kokkos_ENABLE_AGGRESSIVE_VECTORIZATION ON)
  else()
    set(Kokkos_ENABLE_DEBUG TRUE)
    set(Kokkos_ENABLE_AGGRESSIVE_VECTORIZATION OFF)
  endif()
  if (HAERO_ENABLE_GPU)
    if (HAERO_AMD_GPU)
      set(Kokkos_ENABLE_HIP ON CACHE BOOL "Enable HIP Kokkos backend")
    elseif (HAERO_INTEL_GPU)
      set(Kokkos_ENABLE_SYCL ON CACHE BOOL "Enable SYCL Kokkos backend")
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
    set(Kokkos_ENABLE_CUDA OFF CACHE BOOL "Enable CUDA Kokkos backend")
    set(Kokkos_ENABLE_HIP OFF CACHE BOOL "Enable HIP Kokkos backend")
    if (OPENMP_FOUND)
      message(STATUS "OpenMP is enabled")
      set(Kokkos_ENABLE_OPENMP ON CACHE BOOL "Enable OpenMP Kokkos backend")
    else()
      message(STATUS "OpenMP is not supported. Using threads backend instead.")
      set(Kokkos_ENABLE_OPENMP OFF CACHE BOOL "Enable OpenMP Kokkos backend")
      set(Kokkos_ENABLE_PTHREAD ON CACHE BOOL "Enable pthreads Kokkos backend")
    endif()
  endif()

  # Set EKAT config options
  set(EKAT_ENABLE_MPI ${HAERO_ENABLE_MPI} CACHE BOOL "Enable MPI")
  set(EKAT_ENABLE_TESTS OFF CACHE BOOL "Whether tests should be built")
  option(EKAT_ENABLE_KOKKOS       "Enable EKAT kokkos utilities" ON)
  option(EKAT_ENABLE_YAML_PARSER  "Enable EKAT yaml parsing utilities" ON)
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/ekat ${CMAKE_BINARY_DIR}/externals/ekat)
endif()

add_library(ext_libraries INTERFACE)
target_link_libraries (ext_libraries INTERFACE ekat::AllLibs;${HAERO_EXT_LIBRARIES})

# Add all the libraries to the list maintained for the standalone build config file.
set(HAERO_LIBRARIES ekat::AllLibs;${HAERO_EXT_LIBRARIES};${HAERO_LIBRARIES} PARENT_SCOPE)
list(REMOVE_DUPLICATES HAERO_EXT_INCLUDE_DIRS)
set(HAERO_EXT_INCLUDE_DIRS ${HAERO_EXT_INCLUDE_DIRS} PARENT_SCOPE)

#----------------------
# Installation targets
#----------------------

install(DIRECTORY ${PROJECT_BINARY_DIR}/include/ DESTINATION include)
install(DIRECTORY ${PROJECT_BINARY_DIR}/bin/ DESTINATION bin)
