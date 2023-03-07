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
  # Copy EKAT's test-launcher and nvcc_wrapper into place.
  file(COPY ${EKAT_SOURCE_DIR}/bin/test-launcher DESTINATION ${CMAKE_BINARY_DIR}/bin)
  if (HAERO_ENABLE_GPU)
    file(COPY ${EKAT_SOURCE_DIR}/bin/nvcc_wrapper DESTINATION ${CMAKE_BINARY_DIR}/bin)
  endif()
  list(APPEND EKAT_INCLUDE_DIRS
       ${EKAT_SOURCE_DIR}/src
       ${EKAT_BINARY_DIR}/src
       ${EKAT_SOURCE_DIR}/extern/kokkos/tpls/desul/include
       ${EKAT_SOURCE_DIR}/extern/kokkos/core/src
       ${EKAT_SOURCE_DIR}/extern/kokkos/containers/src
       ${EKAT_BINARY_DIR}/../kokkos
       ${EKAT_SOURCE_DIR}/extern/spdlog/include
       ${EKAT_SOURCE_DIR}/extern/yaml-cpp/include)
  set(EKAT_TPL_LIBRARIES kokkos;yaml-cpp;spdlog)
  # FIXME: this is gross, and should be replaced when EKAT has a mechanism for
  # FIXME: describing its installed configuration
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

# FIXME: The following is broken and must either be fixed or removed.
# FIXME: Specifically, it assumes that EKAT can be installed as an external
# FIXME: package instead of being built in-tree. At the time of writing
# FIXME: (1/4/2023), I don't believe we have figured out how to build ekat
# FIXME: as an external package in a way that properly sets up compilers and
# FIXME: flags, because of Kokkos hokkos pokkos. -JNJ
if (HAERO_ENABLE_CHEMISTRY)
  message(FATAL_ERROR "HAERO_ENABLE_CHEMISTRY is unsupported at the moment!")

  #----------------------------------------------------------------------------
  # OpenBLAS linear algebra library: needed by Tines on CPU builds
  # (We can skip it if we use Intel's compilers, which provide MKL.)
  #----------------------------------------------------------------------------
  if (NOT HAERO_ENABLE_GPU AND NOT CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    add_library(openblas STATIC IMPORTED GLOBAL)
    set_target_properties(openblas PROPERTIES IMPORTED_LOCATION ${OPENBLAS_LIBRARY_DIR}/libopenblas.a)
    if (NOT EXISTS ${OPENBLAS_LIBRARY})
      set(OPENBLAS_OPTS libs netlib shared
                        CC=${CMAKE_C_COMPILER}
                        FC=${CMAKE_Fortran_COMPILER})
      if (NOT APPLE)
        set(OPENBLAS_OPTS ${OPENBLAS_OPTS} -j)
      endif()
      ExternalProject_Add(openblas_proj
                          PREFIX ${CMAKE_CURRENT_BINARY_DIR}/OpenBLAS
                          SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/OpenBLAS
                          BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/OpenBLAS
                          INSTALL_DIR ${PROJECT_BINARY_DIR}
                          CONFIGURE_COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/OpenBLAS ${CMAKE_CURRENT_BINARY_DIR}/OpenBLAS
                          BUILD_COMMAND ${MAKE} ${OPENBLAS_OPTS}
                          LOG_BUILD TRUE
                          INSTALL_COMMAND ${MAKE} PREFIX=${PROJECT_BINARY_DIR} install
                          LOG_INSTALL TRUE)

      # Because the OpenBLAS people feel that life isn't complicated enough, they
      # insist on creating libopenblas.a as a symlink. Here we crudely undo this
      # needlessly overcomplicated operation.
      ExternalProject_Add_Step(openblas_proj remove_stupid_openblas_symlink
        COMMAND ${CMAKE_COMMAND} -E rename libopenblas.a libopenblas.a.symlink
        COMMAND ${CMAKE_COMMAND} -E copy libopenblas.a.symlink libopenblas.a
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/lib
        DEPENDEES install)

      add_dependencies(openblas openblas_proj)
    endif()
    list(APPEND HAERO_EXT_INCLUDE_DIRS ${OPENBLAS_INCLUDE_DIR})
    if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
      set(HAERO_EXT_LIBRARIES ${HAERO_EXT_LIBRARIES};gfortran)
    endif()
    set(HAERO_EXT_LIBRARIES openblas;${HAERO_EXT_LIBRARIES})
  endif()

  #----------------------------------------------------------------------------
  # Tines solver library
  #----------------------------------------------------------------------------
  if (TINES_LIBRARY MATCHES ".a")
    add_library(tines STATIC IMPORTED GLOBAL)
  else()
    add_library(tines SHARED IMPORTED GLOBAL)
  endif()
  set_target_properties(tines PROPERTIES IMPORTED_LOCATION ${TINES_LIBRARY})
  if (NOT EXISTS ${TINES_LIBRARY})
    if (NOT HAERO_ENABLE_GPU AND NOT CMAKE_CXX_COMPILER_ID MATCHES "Intel")
      set(TINES_DEPENDS openblas ekat)
    else()
      set(TINES_DEPENDS ekat)
    endif()
    string(REPLACE "${CMAKE_INSTALL_LIBDIR}" "" KOKKOS_INSTALL_PATH ${EKAT_LIBRARY_DIR})
    set(TINES_CXX_FLAGS "-I${EKAT_LIBRARY_DIR}/../include/kokkos -I${PROJECT_BINARY_DIR}/include")
    set(TINES_CMAKE_OPTS -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}
                         -DCMAKE_INSTALL_LIBDIR=${CMAKE_INSTALL_LIBDIR}
                         -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                         -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                         -DCMAKE_CXX_FLAGS=${TINES_CXX_FLAGS}
                         -DTINES_ENABLE_VERBOSE=OFF
                         -DTINES_ENABLE_TEST=OFF
                         -DTINES_ENABLE_EXAMPLE=OFF
                         -DKOKKOS_INSTALL_PATH=${KOKKOS_INSTALL_PATH}
                         -DYAML_INSTALL_PATH=${EKAT_INSTALL_ROOT}/share
                       )
    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
      set(TINES_CMAKE_OPTS ${TINES_CMAKE_OPTS}
                           -DTINES_ENABLE_DEBUG=ON)
    else()
      set(TINES_CMAKE_OPTS ${TINES_CMAKE_OPTS}
                           -DTINES_ENABLE_DEBUG=OFF)
    endif()
    if (HAERO_ENABLE_GPU)
      # We use EKAT's C++ compiler helper script for CUDA.
      set(TINES_CMAKE_OPTS ${TINES_CMAKE_OPTS}
          -DCMAKE_CXX_COMPILER=${EKAT_NVCC_WRAPPER})
    else()
      set(TINES_CMAKE_OPTS ${TINES_CMAKE_OPTS}
                           -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER})
      if (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
        set(TINES_CMAKE_OPTS ${TINES_CMAKE_OPTS} -DTINES_ENABLE_MKL=ON)
      else()
        # We need OpenBLAS for non-GPU builds that don't use Intel's compiler.
        # Extract path to openblas install directory from OPENBLAS_LIBRARY_DIR
        # variable for Tines.
        set(OPENBLAS_INSTALL_PATH ${OPENBLAS_LIBRARY_DIR}/..)
        set(TINES_CMAKE_OPTS ${TINES_CMAKE_OPTS}
                             -DOPENBLAS_INSTALL_PATH=${OPENBLAS_INSTALL_PATH})
      endif()
    endif()

    ExternalProject_Add(tines_proj
                        PREFIX ${CMAKE_CURRENT_BINARY_DIR}/Tines
                        SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/Tines/src
                        BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/Tines
                        INSTALL_DIR ${PROJECT_BINARY_DIR}
                        CMAKE_ARGS ${TINES_CMAKE_OPTS}
                        DEPENDS ${TINES_DEPENDS}
                        LOG_CONFIGURE TRUE
                        BUILD_COMMAND ${MAKE} -j
                        LOG_BUILD TRUE
                        INSTALL_COMMAND ${MAKE} install
                        LOG_INSTALL TRUE)

    # Apparently Tines also leaves a config file that confuses TChem, so we
    # delete it.
    ExternalProject_Add_Step(tines_proj delete_bad_tines_config
      COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_CURRENT_BINARY_DIR}/Tines/core/TinesConfig.cmake
      DEPENDEES install)

    add_dependencies(tines tines_proj)
  endif()
  list(APPEND HAERO_EXT_INCLUDE_DIRS ${TINES_INCLUDE_DIR})
  set(HAERO_EXT_LIBRARIES tines;${HAERO_EXT_LIBRARIES})

  #----------------------------------------------------------------------------
  # TChem chemical kinetics library
  #----------------------------------------------------------------------------
  if (TCHEM_LIBRARY MATCHES ".a")
    add_library(tchem STATIC IMPORTED GLOBAL)
  else()
    add_library(tchem SHARED IMPORTED GLOBAL)
  endif()
  set_target_properties(tchem PROPERTIES IMPORTED_LOCATION ${TCHEM_LIBRARY})
  set(TCHEM_CXX_FLAGS "-I${EKAT_LIBRARY_DIR}/../include/kokkos -I${PROJECT_BINARY_DIR}/include")
  if (NOT EXISTS ${TCHEM_LIBRARY})
    string(REPLACE "${CMAKE_INSTALL_LIBDIR}" "" KOKKOS_INSTALL_PATH ${EKAT_LIBRARY_DIR})
    string(REPLACE "${CMAKE_INSTALL_LIBDIR}" "" TINES_INSTALL_PATH ${TINES_LIBRARY_DIR})
    set(TCHEM_CMAKE_OPTS -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}
                         -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                         -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                         -DTCHEM_INSTALL_LIB_PATH=${CMAKE_INSTALL_LIBDIR}
                         -DTCHEM_ENABLE_VERBOSE=OFF
                         -DTCHEM_ENABLE_KOKKOS=ON
                         -DTCHEM_ENABLE_TEST=OFF
                         -DTCHEM_ENABLE_EXAMPLE=ON
                         -DKOKKOS_INSTALL_PATH=${KOKKOS_INSTALL_PATH}
                         -DTINES_INSTALL_PATH=${TINES_INSTALL_PATH}
                         -DTCHEM_ENABLE_TPL_YAML_CPP=ON
                         # Note(mjs): this enables automatic differentiation, and
                         # is currently only used by the chem_driver unit test.
                         # Keeping it because it seems to be low-impact to compile
                         # time and makes the test run fast
                         -DTCHEM_ENABLE_SACADO_JACOBIAN_ATMOSPHERIC_CHEMISTRY=ON
                       )
    if (HAERO_ENABLE_GPU)
      # We use EKAT's C++ compiler helper script. Also, we need to nudge TChem
      # toward the kokkos/ subdirectory of our include path in this case (for some
      # reason).
      set(TCHEM_CMAKE_OPTS ${TCHEM_CMAKE_OPTS}
          -DCMAKE_CXX_COMPILER=${EKAT_NVCC_WRAPPER})
    else()
      set(TCHEM_CMAKE_OPTS ${TCHEM_CMAKE_OPTS}
          -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER})
    endif()

    set(TCHEM_CMAKE_OPTS ${TCHEM_CMAKE_OPTS} -DCMAKE_CXX_FLAGS=${TCHEM_CXX_FLAGS})
    ExternalProject_Add(tchem_proj
                        PREFIX ${CMAKE_CURRENT_BINARY_DIR}/TChem
                        SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/TChem/src
                        BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/TChem
                        INSTALL_DIR ${PROJECT_BINARY_DIR}
                        CMAKE_ARGS ${TCHEM_CMAKE_OPTS}
                        DEPENDS tines ekat
                        LOG_CONFIGURE TRUE
                        BUILD_COMMAND ${MAKE} -j
                        LOG_BUILD TRUE
                        INSTALL_COMMAND ${MAKE} install
                        LOG_INSTALL TRUE)
    add_dependencies(tchem tchem_proj)
  endif()
  list(APPEND HAERO_EXT_INCLUDE_DIRS ${TCHEM_INCLUDE_DIR})
  set(HAERO_EXT_LIBRARIES tchem;${HAERO_EXT_LIBRARIES})
endif()

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

if (HAERO_ENABLE_CHEMISTRY)
  if (NOT HAERO_ENABLE_GPU AND OPENBLAS_LIBRARY_DIR)
    if (${OPENBLAS_LIBRARY_DIR} MATCHES ${PROJECT_BINARY_DIR}) # we built OpenBLAS
      install(FILES ${OPENBLAS_LIBRARY} DESTINATION ${CMAKE_INSTALL_LIBDIR})
    endif()
  endif()
  if (${TINES_LIBRARY_DIR} MATCHES ${PROJECT_BINARY_DIR}) # we built Tines
    install(FILES ${TINES_LIBRARY} DESTINATION ${CMAKE_INSTALL_LIBDIR})
  endif()
  if (${TCHEM_LIBRARY_DIR} MATCHES ${PROJECT_BINARY_DIR}) # we built tchem
    install(FILES ${TCHEM_LIBRARY} DESTINATION ${CMAKE_INSTALL_LIBDIR})
  endif()
endif()
