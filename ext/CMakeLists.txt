include(ExternalProject)
include_directories(${PROJECT_BINARY_DIR}/include)

# Create canonical top-level directories in the build directory.
foreach(dir include;lib;bin;share)
  set(top_level_dir ${PROJECT_BINARY_DIR}/${dir})
  if (NOT EXISTS ${top_level_dir})
    execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${top_level_dir})
  endif()
endforeach()

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

    # Build zlib.
    set(ZLIB_CONFIG_OPTS --prefix=${PROJECT_BINARY_DIR} --static)
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

    set(HDF5_ENV CC=${CMAKE_C_COMPILER})
    set(HDF5_CONFIG_OPTS --prefix=${PROJECT_BINARY_DIR} --disable-shared
        --with-zlib=${PROJECT_BINARY_DIR}/include,${PROJECT_BINARY_DIR}/lib
        --disable-fortran --disable-cxx --disable-java --disable-tools
        --disable-tests)
    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
      set(HDF5_CONFIG_OPTS ${HDF5_CONFIG_OPTS} --enable-build-mode=debug)
    else()
      set(HDF5_CONFIG_OPTS ${HDF5_CONFIG_OPTS} --enable-build-mode=production)
    endif()
    if (HAVE_MPI)
      set(HDF5_CONFIG_OPTS ${HDF5_CONFIG_OPTS} --enable-parallel)
    endif()
    ExternalProject_Add(hdf5_proj
                        PREFIX ${CMAKE_CURRENT_BINARY_DIR}/hdf5
                        SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/hdf5
                        BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/hdf5
                        INSTALL_DIR ${PROJECT_BINARY_DIR}
                        CONFIGURE_COMMAND env CC=${CMAKE_C_COMPILER} CFLAGS=${CMAKE_C_FLAGS} ${CMAKE_CURRENT_SOURCE_DIR}/hdf5/configure ${HDF5_CONFIG_OPTS}
                        DEPENDS zlib_proj
                        LOG_CONFIGURE TRUE
                        BUILD_COMMAND make -j
                        LOG_BUILD TRUE
                        INSTALL_COMMAND make install
                        LOG_INSTALL TRUE)
    add_dependencies(hdf5 hdf5_proj)
  else()
    message(FATAL_ERROR "Couldn't find HDF5 library: ${HDF5_LIBRARY}")
  endif()
else()
  if (NOT ${HDF5_LIBRARY} MATCHES ${PROJECT_BINARY_DIR}) # if we didn't build this ourselves...
    # Find libz.
    find_package(ZLIB REQUIRED)
    set(HAERO_TP_LIBS ${ZLIB_LIBRARIES};${HAERO_TP_LIBS})
  endif()
endif()
set(HAERO_TP_LIBS hdf5_hl;hdf5;dl;z;${HAERO_TP_LIBS})
set(HAERO_TP_INCDIRS ${HAERO_TP_INCDIRS};${HDF5_INCLUDE_DIR})

#----------------------------------------------------------------------------
# NetCDF C library.
#----------------------------------------------------------------------------
add_library(netcdf STATIC IMPORTED GLOBAL)
set_target_properties(netcdf PROPERTIES IMPORTED_LOCATION ${NETCDF_LIBRARY})
if (NOT EXISTS ${NETCDF_LIBRARY})
  set(NETCDF_CONFIG_OPTS --prefix=${PROJECT_BINARY_DIR} --disable-shared
      --enable-netcdf4 --disable-dap --disable-testsets --disable-examples)
  # HDF5 headers should be given in CPPFLAGS, but the NetCDF folks aren't very
  # careful about this stuff.
  set(NETCDF_C_FLAGS "${CMAKE_C_FLAGS} -I${HDF5_INCLUDE_DIR}")
  ExternalProject_Add(netcdf_proj
                      PREFIX ${CMAKE_CURRENT_BINARY_DIR}/netcdf-c
                      SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/netcdf-c
                      BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/netcdf-c
                      INSTALL_DIR ${PROJECT_BINARY_DIR}
                      CMAKE_ARGS ${NETCDF_CMAKE_OPTS}
                      DEPENDS hdf5
                      CONFIGURE_COMMAND env CC=${CMAKE_C_COMPILER} CFLAGS=${NETCDF_C_FLAGS} LDFLAGS=-L${HDF5_LIBRARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/netcdf-c/configure ${NETCDF_CONFIG_OPTS}
                      LOG_CONFIGURE TRUE
                      BUILD_COMMAND make -j
                      LOG_BUILD TRUE
                      INSTALL_COMMAND make install
                      LOG_INSTALL TRUE)
  add_dependencies(netcdf netcdf_proj)
endif()
set(HAERO_TP_LIBS netcdf;${HAERO_TP_LIBS})
set(HAERO_TP_INCDIRS ${HAERO_TP_INCDIRS};${NETCDF_INCLUDE_DIR})

# NetCDF Fortran library.
add_library(netcdff STATIC IMPORTED GLOBAL)
set_target_properties(netcdff PROPERTIES IMPORTED_LOCATION ${NETCDFF_LIBRARY})
if (NOT EXISTS ${NETCDFF_LIBRARY})
  set(NETCDFF_CONFIG_OPTS --prefix=${PROJECT_BINARY_DIR})
  get_filename_component(NETCDFF_LIBHDF5 ${HDF5_LIBRARY} NAME_WLE)
  string(REPLACE "lib" "" NETCDFF_LIBHDF5 ${NETCDFF_LIBHDF5})
  get_filename_component(NETCDFF_LIBHDF5_HL ${HDF5_HL_LIBRARY} NAME_WLE)
  string(REPLACE "lib" "" NETCDFF_LIBHDF5_HL ${NETCDFF_LIBHDF5_HL})
  set(NETCDFF_LIBS "-l${NETCDFF_LIBHDF5_HL} -l${NETCDFF_LIBHDF5} -lz -ldl")
  ExternalProject_Add(netcdff_proj
                      PREFIX ${CMAKE_CURRENT_BINARY_DIR}/netcdf-fortran
                      SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/netcdf-fortran
                      BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/netcdf-fortran
                      INSTALL_DIR ${PROJECT_BINARY_DIR}
                      CMAKE_ARGS ${NETCDFF_CMAKE_OPTS}
                      DEPENDS netcdf
                      CONFIGURE_COMMAND env CC=${CMAKE_C_COMPILER} FC=${CMAKE_Fortran_COMPILER} CFLAGS=${CMAKE_C_FLAGS} FFLAGS=${CMAKE_Fortran_FLAGS} CPPFLAGS=-I${NETCDF_INCLUDE_DIR} LDFLAGS=-L${NETCDF_LIBRARY_DIR} LIBS=${NETCDFF_LIBS} ${CMAKE_CURRENT_SOURCE_DIR}/netcdf-fortran/configure ${NETCDFF_CONFIG_OPTS}
                      LOG_CONFIGURE TRUE
                      BUILD_COMMAND make -j
                      LOG_BUILD TRUE
                      INSTALL_COMMAND make install
                      LOG_INSTALL TRUE)
  add_dependencies(netcdff netcdff_proj)
endif()
set(HAERO_TP_LIBS netcdff;${HAERO_TP_LIBS})
set(HAERO_TP_INCDIRS ${HAERO_TP_INCDIRS};${NETCDFF_INCLUDE_DIR})

#----------------------------------------------------------------------------
# yaml-cpp C++ YAML library.
#----------------------------------------------------------------------------
if (YAMLCPP_LIBRARY MATCHES ".a")
  add_library(yaml_cpp STATIC IMPORTED GLOBAL)
else()
  add_library(yaml_cpp SHARED IMPORTED GLOBAL)
endif()
set_target_properties(yaml_cpp PROPERTIES IMPORTED_LOCATION ${YAMLCPP_LIBRARY_DIR}/${YAMLCPP_LIBRARY})
if (YAMLCPP_LIBRARY_DIR STREQUAL ${PROJECT_BINARY_DIR}/lib)
  ExternalProject_Add(yaml_cpp_proj
                      PREFIX ${CMAKE_CURRENT_BINARY_DIR}/yaml-cpp
                      SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/yaml-cpp
                      BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/yaml-cpp
                      INSTALL_DIR ${PROJECT_BINARY_DIR}
                      CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}
                      LOG_CONFIGURE TRUE
                      BUILD_COMMAND make -j
                      LOG_BUILD TRUE
                      INSTALL_COMMAND make install
                      LOG_INSTALL TRUE)
  add_dependencies(yaml_cpp yaml_cpp_proj)
endif()
set(HAERO_TP_LIBS yaml_cpp;${HAERO_TP_LIBS})
list(FIND HAERO_TP_INCDIRS ${YAMLCPP_INCLUDE_DIR} YAMLCPP_INDEX)
if (YAMLCPP_INDEX EQUAL -1)
  set(HAERO_TP_INCDIRS ${HAERO_TP_INCDIRS};${YAMLCPP_INCLUDE_DIR})
endif()

# Add all the libraries to the build system at large.
set(HAERO_LIBRARIES ${HAERO_TP_LIBS};${HAERO_LIBRARIES} PARENT_SCOPE)
set(HAERO_INCDIRS ${HAERO_INCDIRS};${HAERO_TP_INCDIRS} PARENT_SCOPE)

# Installation targets
if (${HDF5_LIBRARY} MATCHES ${PROJECT_BINARY_DIR}) # we built HDF5
  add_library(z STATIC IMPORTED GLOBAL)
  set_target_properties(z PROPERTIES IMPORTED_LOCATION ${PROJECT_BINARY_DIR}/lib/libz.a)
  install(FILES ${PROJECT_BINARY_DIR}/lib/libz.a ${HDF5_LIBRARY} ${HDF5_HL_LIBRARY} DESTINATION lib)
endif()
if (${NETCDF_LIBRARY} MATCHES ${PROJECT_BINARY_DIR}) # we built NetCDF-C
  install(FILES ${NETCDF_LIBRARY} DESTINATION lib)
endif()
if (${NETCDFF_LIBRARY} MATCHES ${PROJECT_BINARY_DIR}) # we built NetCDF-Fortran
  install(FILES ${NETCDFF_LIBRARY} DESTINATION lib)
endif()
