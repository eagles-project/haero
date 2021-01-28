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
  message(STATUS "Updating git submodules recursively...")
  if (NOT SUBMODULES_UPDATED)
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
find_library(Z_LIBRARY libz.a)
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
                      BUILD_COMMAND make -j
                      LOG_BUILD TRUE
                      INSTALL_COMMAND make install
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
    message(STATUS "HDF5 not found on system. Building locally.")
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
set(HAERO_EXT_DRIVER_LIBRARIES hdf5_hl;hdf5;z;${HAERO_EXT_DRIVER_LIBRARIES})
set(HAERO_EXT_DRIVER_INCDIRS ${HAERO_EXT_DRIVER_INCDIRS};${HDF5_INCLUDE_DIR})

#----------------------------------------------------------------------------
# NetCDF C library.
#----------------------------------------------------------------------------
add_library(netcdf STATIC IMPORTED GLOBAL)
set_target_properties(netcdf PROPERTIES IMPORTED_LOCATION ${NETCDF_LIBRARY})
if (NOT EXISTS ${NETCDF_LIBRARY})
  message(STATUS "NetCDF not found on system. Building locally.")
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
                      BUILD_COMMAND make
                      LOG_INSTALL TRUE)
  add_dependencies(netcdf netcdf_proj)
endif()
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
set_target_properties(ekat PROPERTIES IMPORTED_LOCATION ${EKAT_LIBRARY_DIR}/${EKAT_LIBRARY})
if (YAMLCPP_LIBRARY MATCHES ".a")
  add_library(yaml_cpp STATIC IMPORTED GLOBAL)
else()
  add_library(yaml_cpp SHARED IMPORTED GLOBAL)
endif()
set_target_properties(yaml_cpp PROPERTIES IMPORTED_LOCATION ${YAMLCPP_LIBRARY_DIR}/${YAMLCPP_LIBRARY})

# Check if we're supposed to build ekat internally
get_filename_component(EKAT_INSTALL_ROOT ${EKAT_LIBRARY_DIR}/.. ABSOLUTE)
if (EKAT_INSTALL_ROOT STREQUAL PROJECT_BINARY_DIR)
  message(STATUS "EKAT library not found on system. Building locally.")
  update_submodules()
  ExternalProject_Add(ekat_proj
                      PREFIX ${CMAKE_CURRENT_BINARY_DIR}/ekat
                      SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ekat
                      BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/ekat
                      INSTALL_DIR ${PROJECT_BINARY_DIR}
                      CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR} -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                      LOG_CONFIGURE TRUE
                      BUILD_COMMAND make -j
                      LOG_BUILD TRUE
                      INSTALL_COMMAND make install
                      LOG_INSTALL TRUE)
  add_dependencies(yaml_cpp ekat_proj)
  add_dependencies(ekat ekat_proj)
endif()
add_library(kokkoscore STATIC IMPORTED GLOBAL)
set_target_properties(kokkoscore PROPERTIES IMPORTED_LOCATION ${EKAT_LIBRARY_DIR}/libkokkoscore.a)
add_library(kokkoscontainers STATIC IMPORTED GLOBAL)
set_target_properties(kokkoscontainers PROPERTIES IMPORTED_LOCATION ${EKAT_LIBRARY_DIR}/libkokkoscontainers.a)
set(HAERO_EXT_LIBRARIES ekat;kokkoscontainers;kokkoscore;${HAERO_EXT_LIBRARIES})
set(HAERO_EXT_DRIVER_LIBRARIES yaml_cpp;${HAERO_EXT_DRIVER_LIBRARIES})
list(FIND HAERO_EXT_INCDIRS ${EKAT_INCLUDE_DIR} EKAT_INDEX)
if (EKAT_INDEX EQUAL -1)
  set(HAERO_EXT_INCDIRS ${HAERO_EXT_INCDIRS};${EKAT_INCLUDE_DIR};${EKAT_INCLUDE_DIR}/ekat;${EKAT_INCLUDE_DIR}/kokkos)
endif()
list(FIND HAERO_EXT_INCDIRS ${YAMLCPP_INCLUDE_DIR} YAMLCPP_INDEX)
if (YAMLCPP_INDEX EQUAL -1)
  set(HAERO_EXT_INCDIRS ${HAERO_EXT_INCDIRS};${YAMLCPP_INCLUDE_DIR})
endif()

#----------------------------------------------------------------------------
# Kokkos-Kernels library
#----------------------------------------------------------------------------
if (KOKKOS_KERNELS_LIBRARY MATCHES ".a")
  add_library(kokkoskernels STATIC IMPORTED GLOBAL)
else()
  add_library(kokkoskernels SHARED IMPORTED GLOBAL)
endif()
set_target_properties(kokkoskernels PROPERTIES IMPORTED_LOCATION ${KOKKOS_KERNELS_LIBRARY_DIR}/${KOKKOS_KERNELS_LIBRARY})

if (NOT EXISTS ${KOKKOS_KERNELS_LIBRARY})
  message(STATUS "kokkos-kernels library not found on system. Building locally.")
  update_submodules()
  string(REPLACE "/lib" "" KOKKOS_DIR ${EKAT_LIBRARY_DIR})
  set(KK_CMAKE_OPTS -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}
                    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                    -DKokkos_ROOT=${KOKKOS_DIR}
                    -DKokosKernels_INST_LAYOUTRIGHT:BOOL=ON
                    -DCMAKE_CXX_EXTENSIONS:BOOL=OFF
                    )
  ExternalProject_Add(kokkoskernels_proj
                      PREFIX ${CMAKE_CURRENT_BINARY_DIR}/kokkos-kernels
                      SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/kokkos-kernels
                      BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/kokkos-kernels
                      INSTALL_DIR ${PROJECT_BINARY_DIR}
                      CMAKE_ARGS ${KK_CMAKE_OPTS}
                      DEPENDS ekat
                      LOG_CONFIGURE TRUE
                      BUILD_COMMAND make -j
                      LOG_BUILD TRUE
                      INSTALL_COMMAND make install
                      LOG_INSTALL TRUE)
  add_dependencies(kokkoskernels kokkoskernels_proj)

  # Our version of Kokkos has an install process that's kind of broken:
  # 1. It leaves a KokkosConfig.cmake file in its build directory in addition to
  #    the installation path, and this confuses both kokkos-kernels and TChem.
  # 2. For some reason, Kokkos_Sort.hpp and Kokkos_Random.hpp are left out of
  #    the kokkos/ include path. This oddity likely indicates a version
  #    mismatch.
  # This custom step attempts to address these issues.
  if (EKAT_INSTALL_ROOT STREQUAL PROJECT_BINARY_DIR)
    ExternalProject_Add_Step(kokkoskernels_proj delete_bad_kokkos_config
      COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_CURRENT_BINARY_DIR}/ekat/externals/kokkos/KokkosConfig.cmake
      DEPENDS ekat
      DEPENDERS configure)
  endif()
endif()
set(HAERO_EXT_LIBRARIES kokkoskernels;${HAERO_EXT_LIBRARIES})

#----------------------------------------------------------------------------
# TChem chemical kinetics library
#----------------------------------------------------------------------------
if (TCHEM_LIBRARY MATCHES ".a")
  add_library(tchem STATIC IMPORTED GLOBAL)
else()
  add_library(tchem SHARED IMPORTED GLOBAL)
endif()
set_target_properties(tchem PROPERTIES IMPORTED_LOCATION ${TCHEM_LIBRARY_DIR}/${TCHEM_LIBRARY})

if (NOT EXISTS ${TCHEM_LIBRARY})
  message(STATUS "TChem library not found on system. Building locally.")
  update_submodules()
  string(REPLACE "/lib" "" KOKKOS_INSTALL_PATH ${EKAT_LIBRARY_DIR})
  set(TCHEM_CMAKE_OPTS -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}
                       -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                       -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                       -DTCHEM_ENABLE_VERBOSE=ON
                       -DTCHEM_ENABLE_KOKKOS=ON
                       -DTCHEM_ENABLE_KOKKOSKERNELS=ON
                       -DTCHEM_ENABLE_TEST=OFF
                       -DTCHEM_ENABLE_EXAMPLE=ON
                       -DKOKKOS_INSTALL_PATH=${KOKKOS_INSTALL_PATH}
                     )
  ExternalProject_Add(tchem_proj
                      PREFIX ${CMAKE_CURRENT_BINARY_DIR}/TChem
                      SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/TChem/src
                      BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/TChem
                      INSTALL_DIR ${PROJECT_BINARY_DIR}
                      CMAKE_ARGS ${TCHEM_CMAKE_OPTS}
                      DEPENDS kokkoskernels
                      LOG_CONFIGURE TRUE
                      BUILD_COMMAND make -j
                      LOG_BUILD TRUE
                      INSTALL_COMMAND make install
                      LOG_INSTALL TRUE)
  add_dependencies(tchem tchem_proj)
endif()
set(HAERO_EXT_LIBRARIES tchem;${HAERO_EXT_LIBRARIES})

# Add all the libraries to the build system at large.
set(HAERO_LIBRARIES ${HAERO_EXT_LIBRARIES};${HAERO_LIBRARIES} PARENT_SCOPE)
set(HAERO_INCDIRS ${HAERO_INCDIRS};${HAERO_EXT_INCDIRS} PARENT_SCOPE)
set(HAERO_DRIVER_LIBRARIES ${HAERO_EXT_DRIVER_LIBRARIES} PARENT_SCOPE)
set(HAERO_DRIVER_INCDIRS ${HAERO_EXT_DRIVER_INCDIRS} PARENT_SCOPE)

# Installation targets
if (${Z_LIBRARY} MATCHES ${PROJECT_BINARY_DIR}) # we built libz
  install(FILES ${Z_LIBRARY} DESTINATION lib)
endif()
if (${HDF5_LIBRARY} MATCHES ${PROJECT_BINARY_DIR}) # we built HDF5
  install(FILES ${HDF5_LIBRARY} ${HDF5_HL_LIBRARY} DESTINATION lib)
endif()
if (${NETCDF_LIBRARY} MATCHES ${PROJECT_BINARY_DIR}) # we built NetCDF-C
  install(FILES ${NETCDF_LIBRARY} DESTINATION lib)
endif()
if (${EKAT_LIBRARY} MATCHES ${PROJECT_BINARY_DIR}) # we built ekat
  install(FILES ${EKAT_LIBRARY} DESTINATION lib)
endif()
if (${YAMLCPP_LIBRARY} MATCHES ${PROJECT_BINARY_DIR}) # we built yaml-cpp
  install(FILES ${YAMLCPP_LIBRARY} DESTINATION lib)
endif()
if (${TCHEM_LIBRARY} MATCHES ${PROJECT_BINARY_DIR}) # we built tchem
  install(FILES ${TCHEM_LIBRARY} DESTINATION lib)
endif()
