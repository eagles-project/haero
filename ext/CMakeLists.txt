include(ExternalProject)
include(GNUInstallDirs)
include(${PROJECT_SOURCE_DIR}/cmake/cmake_helpers.cmake)
include_directories(${PROJECT_BINARY_DIR}/include)

# Create canonical top-level directories in the build directory.
foreach(dir include;lib;bin;share)
  set(top_level_dir ${PROJECT_BINARY_DIR}/${dir})
  if (NOT EXISTS ${top_level_dir})
    execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${top_level_dir})
  endif()
endforeach()

# Keep track of whether we've updated our submodules.
set(SUBMODULES_UPDATED FALSE)
macro(update_submodules)
  if (NOT SUBMODULES_UPDATED OR HAERO_DISABLE_SUBMODULE_CHECKS)
    message(STATUS "Updating git submodules recursively...")
    execute_process(COMMAND ${GIT} submodule update --init --recursive
                    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
    set(SUBMODULES_UPDATED TRUE)
  endif()
endmacro()

# Many of our dependencies are only needed in standalone mode.
if (HAERO_ENABLE_DRIVER)

#----------------------------------------------------------------------------
# Dynamic loading library
#----------------------------------------------------------------------------
if (NOT APPLE)
  find_library(DL_LIBRARY libdl${CMAKE_SHARED_LIBRARY_SUFFIX} REQUIRED)
  add_library(dl SHARED IMPORTED GLOBAL)
  set_target_properties(dl PROPERTIES IMPORTED_LOCATION ${DL_LIBRARY})
  set(HAERO_EXT_LIBRARIES dl) # Required by kokkoscore.
endif()

#----------------------------------------------------------------------------
# libz compression library
#----------------------------------------------------------------------------
find_library(Z_LIBRARY libz.a libz${CMAKE_SHARED_LIBRARY_SUFFIX})
if (Z_LIBRARY STREQUAL "Z_LIBRARY-NOTFOUND")
  # Build zlib.
  message(STATUS "libz not found on system. Building locally.")
  update_submodules()
  set(ZLIB_CONFIG_OPTS
      --prefix=${PROJECT_BINARY_DIR}
      --static
      --libdir=${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
  ExternalProject_Add(zlib_proj
                      PREFIX ${CMAKE_CURRENT_BINARY_DIR}/zlib
                      SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/zlib
                      BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/zlib
                      INSTALL_DIR ${PROJECT_BINARY_DIR}
                      CONFIGURE_COMMAND env CC=${CMAKE_C_COMPILER} CFLAGS=${CMAKE_C_FLAGS} ${CMAKE_CURRENT_SOURCE_DIR}/zlib/configure ${ZLIB_CONFIG_OPTS}
                      LOG_CONFIGURE TRUE
                      BUILD_COMMAND ${MAKE} -j
                      LOG_BUILD TRUE
                      INSTALL_COMMAND ${MAKE} install
                      LOG_INSTALL TRUE)
  set(Z_LIBRARY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}/libz.a)
  set(HDF5_DEPS zlib_proj)
else()
  message(STATUS "Found libz: ${Z_LIBRARY}")
endif()

add_library(z STATIC IMPORTED GLOBAL)
set_target_properties(z PROPERTIES IMPORTED_LOCATION ${Z_LIBRARY})

#----------------------------------------------------------------------------
# HDF5 parallel I/O library
#----------------------------------------------------------------------------
if (HDF5_LIBRARY MATCHES ".so")
  add_library(hdf5 SHARED IMPORTED GLOBAL)
  add_library(hdf5_hl SHARED IMPORTED GLOBAL)
else()
  add_library(hdf5 STATIC IMPORTED GLOBAL)
  add_library(hdf5_hl STATIC IMPORTED GLOBAL)
endif()
set_target_properties(hdf5 PROPERTIES IMPORTED_LOCATION ${HDF5_LIBRARY})
set_target_properties(hdf5_hl PROPERTIES IMPORTED_LOCATION ${HDF5_HL_LIBRARY})
if (NOT EXISTS ${HDF5_LIBRARY})
  if (${HDF5_LIBRARY} MATCHES ${PROJECT_BINARY_DIR})
    update_submodules()

    set(HDF5_CMAKE_OPTS -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}
                        -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                        -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
                        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                        -DZLIB_INCLUDE_DIRS=${PROJECT_BINARY_DIR}/include
                        -DZLIB_LIBARIES=${Z_LIBRARY} -DZLIB_USE_EXTERNAL=OFF
                        -DHDF5_INSTALL_LIB_DIR=${CMAKE_INSTALL_LIBDIR}
                        -DHDF5_ENABLE_Z_LIB_SUPPORT=ON
                        -DHDF5_ENABLE_PARALLEL=${HAVE_MPI}
                        -DHDF5_BUILD_CPP_LIB=OFF
                        -DHDF5_BUILD_FORTRAN=OFF
                        -DHDF5_BUILD_JAVA=OFF
                        -DHDF5_BUILD_TOOLS=OFF
                        -DHDF5_BUILD_UTILS=OFF
                        -DBUILD_TESTING=OFF
                        -DHDF5_BUILD_EXAMPLES=OFF)
    ExternalProject_Add(hdf5_proj
                        PREFIX ${CMAKE_CURRENT_BINARY_DIR}/hdf5
                        SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/hdf5
                        BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/hdf5
                        INSTALL_DIR ${PROJECT_BINARY_DIR}
                        CMAKE_ARGS ${HDF5_CMAKE_OPTS}
                        DEPENDS ${HDF5_DEPS}
                        LOG_CONFIGURE TRUE
                        LOG_BUILD TRUE
                        LOG_INSTALL TRUE)
    add_dependencies(hdf5 hdf5_proj)
  else()
    message(FATAL_ERROR "Couldn't find HDF5 library: ${HDF5_LIBRARY}")
  endif()
else()
  if (NOT ${HDF5_LIBRARY} MATCHES ${PROJECT_BINARY_DIR}) # if we didn't build this ourselves...
    # Find libz.
    find_package(ZLIB REQUIRED)
    set(HAERO_EXT_DRIVER_LIBRARIES ${ZLIB_LIBRARIES};${HAERO_EXT_DRIVER_LIBRARIES})
  endif()
endif()
list(APPEND HAERO_EXT_INCDIRS ${HDF5_INCLUDE_DIR})
set(HAERO_EXT_DRIVER_LIBRARIES hdf5_hl;hdf5;z;${HAERO_EXT_DRIVER_LIBRARIES})
set(HAERO_EXT_DRIVER_INCDIRS ${HAERO_EXT_DRIVER_INCDIRS};${HDF5_INCLUDE_DIR})

#----------------------------------------------------------------------------
# NetCDF C library.
#----------------------------------------------------------------------------
add_library(netcdf STATIC IMPORTED GLOBAL)
set_target_properties(netcdf PROPERTIES IMPORTED_LOCATION ${NETCDF_LIBRARY})
if (NOT EXISTS ${NETCDF_LIBRARY})
  update_submodules()
  set(NETCDF_CMAKE_OPTS -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                        -DCMAKE_C_FLAGS=${NETCDF_C_FLAGS}
                        -DUSE_HDF5=ON
                        -DHDF5_C_LIBRARY=${HDF5_LIBRARY}
                        -DHDF5_HL_LIBRARY=${HDF5_HL_LIBRARY}
                        -DHDF5_C_LIBRARIES=${HDF5_LIBRARY}
                        -DHDF5_HL_LIBRARIES=${HDF5_HL_LIBRARY}
                        -DHDF5_INCLUDE_DIR=${HDF5_INCLUDE_DIR}
                        -DZLIB_LIBRARY=${Z_LIBRARY}
                        -DHAVE_LIBDL=ON
                        -DLIBDL=${DL_LIBRARY}
                        -DENABLE_PARALLEL=${HAVE_MPI}
                        -DHDF5_IS_PARALLEL=${HAVE_MPI}
                        -DHDF5_IS_PARALLEL_MPIO=${HAVE_MPI}
                        -DNC_HAVE_PARALLEL_HDF5=${HAVE_MPI}
                        -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}
                        -DCMAKE_INSTALL_LIBDIR=${CMAKE_INSTALL_LIBDIR}
                        -DENABLE_DAP=OFF
                        -DBUILD_SHARED_LIBS=OFF
                        -DFIND_SHARED_LIBS=OFF
                        -DENABLE_DYNAMIC_LOADING=ON
                        -DENABLE_TESTS=OFF
                        -DBUILD_UTILITIES=ON
                        -DBUILD_EXAMPLES=OFF
                        -DCMAKE_FIND_LIBRARY_SUFFIXES=.a)
  ExternalProject_Add(netcdf_proj
                      PREFIX ${CMAKE_CURRENT_BINARY_DIR}/netcdf-c
                      SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/netcdf-c
                      BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/netcdf-c
                      INSTALL_DIR ${PROJECT_BINARY_DIR}
                      CMAKE_ARGS ${NETCDF_CMAKE_OPTS}
                      DEPENDS hdf5
                      LOG_CONFIGURE TRUE
                      LOG_BUILD TRUE
                      BUILD_COMMAND ${MAKE}
                      LOG_INSTALL TRUE)
  add_dependencies(netcdf netcdf_proj)
endif()
list(APPEND HAERO_EXT_INCDIRS ${NETCDF_INCLUDE_DIR})
set(HAERO_EXT_DRIVER_LIBRARIES netcdf;${HAERO_EXT_DRIVER_LIBRARIES})
set(HAERO_EXT_DRIVER_INCDIRS ${HAERO_EXT_DRIVER_INCDIRS};${NETCDF_INCLUDE_DIR})

endif() # ENABLE_HAERO_DRIVER

#----------------------------------------------------------------------------
# E3SM Kokkos Application Toolkit (EKAT) library and friends.
#----------------------------------------------------------------------------
if (EKAT_LIBRARY MATCHES ".a")
  add_library(ekat STATIC IMPORTED GLOBAL)
else()
  add_library(ekat SHARED IMPORTED GLOBAL)
endif()
set_target_properties(ekat PROPERTIES IMPORTED_LOCATION ${EKAT_LIBRARY})
if (YAMLCPP_LIBRARY MATCHES ".a")
  add_library(yaml_cpp STATIC IMPORTED GLOBAL)
else()
  add_library(yaml_cpp SHARED IMPORTED GLOBAL)
endif()
set_target_properties(yaml_cpp PROPERTIES IMPORTED_LOCATION ${YAMLCPP_LIBRARY})

# Check if we're supposed to build ekat internally
get_filename_component(EKAT_INSTALL_ROOT ${EKAT_LIBRARY_DIR}/.. ABSOLUTE)
if (NOT EXISTS ${EKAT_LIBRARY})
  update_submodules()
  set(EKAT_CMAKE_OPTS -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}
                      -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                      -DCMAKE_INSTALL_LIBDIR=${CMAKE_INSTALL_LIBDIR}
                      -DKokkos_ENABLE_DEPRECATED_CODE=OFF
                      -DKokkos_ENABLE_PROFILING=OFF)
  if (CMAKE_BUILD_TYPE STREQUAL Release)
    set(EKAT_CMAKE_OPTS ${EKAT_CMAKE_OPTS}
                        -DKokkos_ENABLE_DEBUG=FALSE
                        -DKokkos_ENABLE_AGGRESSIVE_VECTORIZATION=ON)
  else()
    set(EKAT_CMAKE_OPTS ${EKAT_CMAKE_OPTS}
                        -DKokkos_ENABLE_DEBUG=TRUE
                        -DKokkos_ENABLE_AGGRESSIVE_VECTORIZATION=OFF)
  endif()
  if (HAERO_DEVICE STREQUAL CPU)
    set(EKAT_CMAKE_OPTS ${EKAT_CMAKE_OPTS}
                        -DKokkos_ENABLE_CUDA=OFF)
    if (NOT APPLE)
      set(EKAT_CMAKE_OPTS ${EKAT_CMAKE_OPTS}
                          -DKokkos_ENABLE_SERIAL=ON)
    endif()
    if (HAVE_OPENMP)
      set(EKAT_CMAKE_OPTS ${EKAT_CMAKE_OPTS}
                          -DKokkos_ENABLE_OPENMP=ON)
    else()
      set(EKAT_CMAKE_OPTS ${EKAT_CMAKE_OPTS}
                          -DKokkos_ENABLE_OPENMP=OFF)
    endif()
  elseif (HAERO_DEVICE STREQUAL CUDA)
    set(KOKKOS_ARCH_FLAG "Kokkos_ARCH_${HAERO_DEVICE_ARCH}")
    set(EKAT_CMAKE_OPTS ${EKAT_CMAKE_OPTS}
                        -DKokkos_ENABLE_SERIAL=ON
                        -DKokkos_ENABLE_CUDA=ON
                        -DKokkos_ENABLE_CUDA_LAMBDA=ON
                        -D${KOKKOS_ARCH_FLAG}:BOOL=ON)
  endif()
  ExternalProject_Add(ekat_proj
                      PREFIX ${CMAKE_CURRENT_BINARY_DIR}/ekat
                      SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ekat
                      BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/ekat
                      INSTALL_DIR ${PROJECT_BINARY_DIR}
                      CMAKE_ARGS ${EKAT_CMAKE_OPTS}
                      LOG_CONFIGURE TRUE
                      BUILD_COMMAND ${MAKE} -j
                      LOG_BUILD TRUE
                      INSTALL_COMMAND ${MAKE} install
                      LOG_INSTALL TRUE)

  # The above install process leaves a KokkosConfig.cmake file in its build
  # directory in addition to the installation path, and this confuses subsequent
  # libraries, so we add a custom step to get rid of it.
  ExternalProject_Add_Step(ekat_proj delete_bad_kokkos_config
    COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_CURRENT_BINARY_DIR}/ekat/externals/kokkos/KokkosConfig.cmake
    DEPENDEES install)

  add_dependencies(yaml_cpp ekat_proj)
  add_dependencies(ekat ekat_proj)
endif()
add_library(kokkoscore STATIC IMPORTED GLOBAL)
set_target_properties(kokkoscore PROPERTIES IMPORTED_LOCATION ${EKAT_LIBRARY_DIR}/libkokkoscore.a)
add_library(kokkoscontainers STATIC IMPORTED GLOBAL)
set_target_properties(kokkoscontainers PROPERTIES IMPORTED_LOCATION ${EKAT_LIBRARY_DIR}/libkokkoscontainers.a)
set(HAERO_EXT_LIBRARIES ekat;kokkoscontainers;kokkoscore;${HAERO_EXT_LIBRARIES})
set(HAERO_EXT_DRIVER_LIBRARIES yaml_cpp;${HAERO_EXT_DRIVER_LIBRARIES})
list(APPEND HAERO_EXT_INCDIRS ${EKAT_INCLUDE_DIR})
list(APPEND HAERO_EXT_INCDIRS ${EKAT_INCLUDE_DIR}/ekat)
list(APPEND HAERO_EXT_INCDIRS ${EKAT_INCLUDE_DIR}/kokkos)
list(APPEND HAERO_EXT_INCDIRS ${YAMLCPP_INCLUDE_DIR})

#----------------------------------------------------------------------------
# OpenBLAS linear algebra library: needed by Tines on CPU builds
# (We can skip it if we use Intel's compilers, which provide MKL.)
#----------------------------------------------------------------------------
if (HAERO_DEVICE STREQUAL "CPU" AND NOT CMAKE_CXX_COMPILER_ID MATCHES "Intel")
  add_library(openblas STATIC IMPORTED GLOBAL)
  set_target_properties(openblas PROPERTIES IMPORTED_LOCATION ${OPENBLAS_LIBRARY_DIR}/libopenblas.a)
  if (NOT EXISTS ${OPENBLAS_LIBRARY})
    update_submodules()
    set(OPENBLAS_OPTS libs CC=${CMAKE_C_COMPILER}
                      FC=${CMAKE_Fortran_COMPILER}
                      NO_SHARED=1 PREFIX=${PROJECT_BINARY_DIR})
    if (NOT APPLE)
      set(OPENBLAS_OPTS ${OPENBLAS_OPTS} -j)
    endif()
    if (HAERO_DOUBLE_PRECISION)
      set(OPENBLAS_OPTS ${OPENBLAS_OPTS} BUILD_DOUBLE=1)
    else()
      set(OPENBLAS_OPTS ${OPENBLAS_OPTS} BUILD_SINGLE=1)
    endif()
    ExternalProject_Add(openblas_proj
                        PREFIX ${CMAKE_CURRENT_BINARY_DIR}/OpenBLAS
                        SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/OpenBLAS
                        BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/OpenBLAS
                        INSTALL_DIR ${PROJECT_BINARY_DIR}
                        CONFIGURE_COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/OpenBLAS ${CMAKE_CURRENT_BINARY_DIR}/OpenBLAS
                        BUILD_COMMAND ${MAKE} ${OPENBLAS_OPTS}
                        LOG_BUILD TRUE
                        INSTALL_COMMAND ${MAKE} ${OPENBLAS_OPTS} install
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
  list(APPEND HAERO_EXT_INCDIRS ${OPENBLAS_INCLUDE_DIR})
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
  update_submodules()
  if (HAERO_DEVICE STREQUAL "CPU" AND NOT CMAKE_CXX_COMPILER_ID MATCHES "Intel")
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
                     )
  if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(TINES_CMAKE_OPTS ${TINES_CMAKE_OPTS}
                         -DTINES_ENABLE_DEBUG=ON)
  else()
    set(TINES_CMAKE_OPTS ${TINES_CMAKE_OPTS}
                         -DTINES_ENABLE_DEBUG=OFF)
  endif()
  if (HAERO_DEVICE STREQUAL "CUDA")
    # We use EKAT's C++ compiler helper script for CUDA.
    set(TINES_CMAKE_OPTS ${TINES_CMAKE_OPTS}
                         -DCMAKE_CXX_COMPILER=${EKAT_LIBRARY_DIR}/../bin/ekat_mpicxx)
  else()
    set(TINES_CMAKE_OPTS ${TINES_CMAKE_OPTS}
                         -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER})
    if (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
      set(TINES_CMAKE_OPTS ${TINES_CMAKE_OPTS}
                           -DTINES_ENABLE_MKL=ON)
    else()
      # We need OpenBLAS for non-GPU builds that don't use Intel's compiler.
      #Extract path to openblas install directory from OPENBLAS_LIBRARY_DIR variable for Tines
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
list(APPEND HAERO_EXT_INCDIRS ${TINES_INCLUDE_DIR})
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
  update_submodules()
  string(REPLACE "${CMAKE_INSTALL_LIBDIR}" "" KOKKOS_INSTALL_PATH ${EKAT_LIBRARY_DIR})
  string(REPLACE "${CMAKE_INSTALL_LIBDIR}" "" TINES_INSTALL_PATH ${TINES_LIBRARY_DIR})
  set(TCHEM_CMAKE_OPTS -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}
                       -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
		       -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                       -DTCHEM_INSTALL_LIB_PATH=${CMAKE_INSTALL_LIBDIR}
                       -DTCHEM_ENABLE_VERBOSE=OFF
                       -DTCHEM_ENABLE_KOKKOS=ON
                       -DTCHEM_ENABLE_TEST=OFF
                       -DTCHEM_ENABLE_EXAMPLE=OFF
                       -DKOKKOS_INSTALL_PATH=${KOKKOS_INSTALL_PATH}
                       -DTINES_INSTALL_PATH=${TINES_INSTALL_PATH}
                     )
  if (HAERO_DEVICE STREQUAL "CUDA")
    # We use EKAT's C++ compiler helper script. Also, we need to nudge TChem
    # toward the kokkos/ subdirectory of our include path in this case (for some
    # reason).
    set(TCHEM_CMAKE_OPTS ${TCHEM_CMAKE_OPTS}
                         -DCMAKE_CXX_COMPILER=${EKAT_LIBRARY_DIR}/../bin/ekat_mpicxx)
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
list(APPEND HAERO_EXT_INCDIRS ${TCHEM_INCLUDE_DIR})
set(HAERO_EXT_LIBRARIES tchem;${HAERO_EXT_LIBRARIES})

# Add all the libraries to the build system at large.
set(HAERO_LIBRARIES ${HAERO_EXT_LIBRARIES};${HAERO_LIBRARIES} PARENT_SCOPE)
list(REMOVE_DUPLICATES HAERO_EXT_INCDIRS)
set(HAERO_INCDIRS ${HAERO_INCDIRS};${HAERO_EXT_INCDIRS} PARENT_SCOPE)
set(HAERO_DRIVER_LIBRARIES ${HAERO_EXT_DRIVER_LIBRARIES} PARENT_SCOPE)
set(HAERO_DRIVER_INCDIRS ${HAERO_EXT_DRIVER_INCDIRS} PARENT_SCOPE)

# Installation targets
install(DIRECTORY ${PROJECT_BINARY_DIR}/include/ DESTINATION include)

if (${Z_LIBRARY_DIR} MATCHES ${PROJECT_BINARY_DIR}) # we built libz
  install(FILES ${Z_LIBRARY} DESTINATION lib)
endif()
if (${HDF5_LIBRARY_DIR} MATCHES ${PROJECT_BINARY_DIR}) # we built HDF5
  install(FILES ${HDF5_LIBRARY} ${HDF5_HL_LIBRARY} DESTINATION lib)
endif()
if (${NETCDF_LIBRARY_DIR} MATCHES ${PROJECT_BINARY_DIR}) # we built NetCDF-C
  install(FILES ${NETCDF_LIBRARY} DESTINATION lib)
endif()
if (${EKAT_LIBRARY_DIR} MATCHES ${PROJECT_BINARY_DIR}) # we built ekat
  install(FILES       ${EKAT_LIBRARY}
                      ${EKAT_LIBRARY_DIR}/libekat_test_main.a
                      ${EKAT_LIBRARY_DIR}/libekat_test_session.a
                      ${EKAT_LIBRARY_DIR}/libkokkoscore.a
                      ${EKAT_LIBRARY_DIR}/libkokkoscontainers.a
          DESTINATION lib)
endif()
if (${YAMLCPP_LIBRARY_DIR} MATCHES ${PROJECT_BINARY_DIR}) # we built yaml-cpp
  install(FILES ${YAMLCPP_LIBRARY} DESTINATION lib)
endif()
if (HAERO_DEVICE STREQUAL "CPU" AND OPENBLAS_LIBRARY_DIR)
  if (${OPENBLAS_LIBRARY_DIR} MATCHES ${PROJECT_BINARY_DIR}) # we built OpenBLAS
    install(FILES ${OPENBLAS_LIBRARY} DESTINATION lib)
  endif()
endif()
if (${TINES_LIBRARY_DIR} MATCHES ${PROJECT_BINARY_DIR}) # we built Tines
  install(FILES ${TINES_LIBRARY} DESTINATION lib)
endif()
if (${TCHEM_LIBRARY_DIR} MATCHES ${PROJECT_BINARY_DIR}) # we built tchem
  install(FILES ${TCHEM_LIBRARY} DESTINATION lib)
endif()
if (${PROJECT_BINARY_DIR}/share)
  install(DIRECTORY ${PROJECT_BINARY_DIR}/share DESTINATION share)
endif()
