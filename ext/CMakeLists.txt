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

# Many of our dependencies are only needed in standalone mode.
if (HAERO_ENABLE_DRIVER)

#----------------------------------------------------------------------------
# Dynamic loading library
#----------------------------------------------------------------------------
find_library(DL_LIBRARY libdl${CMAKE_SHARED_LIBRARY_SUFFIX} REQUIRED)
add_library(dl SHARED IMPORTED GLOBAL)
set_target_properties(dl PROPERTIES IMPORTED_LOCATION ${DL_LIBRARY})

#----------------------------------------------------------------------------
# libz compression library
#----------------------------------------------------------------------------
find_library(Z_LIBRARY libz.a)
if (Z_LIBRARY STREQUAL "Z_LIBRARY-NOTFOUND")
  # Build zlib.
  message(STATUS "libz not found on system. Building locally.")
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
    # We need to build HDF5 locally. First, update git submodules
    execute_process(COMMAND ${GIT} submodule update --init --recursive
                    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})

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
set(HAERO_EXT_DRIVER_LIBRARIES hdf5_hl;hdf5;dl;z;${HAERO_EXT_DRIVER_LIBRARIES})
set(HAERO_EXT_DRIVER_INCDIRS ${HAERO_EXT_DRIVER_INCDIRS};${HDF5_INCLUDE_DIR})

#----------------------------------------------------------------------------
# NetCDF C library.
#----------------------------------------------------------------------------
add_library(netcdf STATIC IMPORTED GLOBAL)
set_target_properties(netcdf PROPERTIES IMPORTED_LOCATION ${NETCDF_LIBRARY})
if (NOT EXISTS ${NETCDF_LIBRARY})
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
