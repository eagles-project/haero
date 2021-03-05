# This macro identifies compilers and third-party library needs
# for particular hosts.
macro(set_up_platform)

  # Are we on Linux?
  if (UNIX AND NOT APPLE)
    set(LINUX ON)
  endif()

  # Do we have bash?
  find_program(BASH bash)
  if (BASH STREQUAL "BASH_NOTFOUND")
    message(FATAL_ERROR "Bash is required, but is not available on this system.")
  endif()

  # Do we have make?
  find_program(MAKE make)
  if (MAKE STREQUAL "MAKE_NOTFOUND")
    message(FATAL_ERROR "Make is required, but is not available on this system.")
  endif()

  # Do we have git?
  find_program(GIT git)
  if (GIT STREQUAL "GIT_NOTFOUND")
    message(WARNING "Git not found. Hope you're not developing on this system.")
    set(HAVE_GIT FALSE)
  else()
    set(HAVE_GIT TRUE)
  endif()

  include(GNUInstallDirs)

  if (HDF5_INCLUDE_DIR)
    if (NOT EXISTS ${HDF5_INCLUDE_DIR})
      message(FATAL_ERROR "Couldn't find HDF5 include dir at ${HDF5_INCLUDE_DIR}.")
    endif()
    message(STATUS "Using HDF5 include dir: ${HDF5_INCLUDE_DIR}.")
  else()
    set(HDF5_INCLUDE_DIR "${CMAKE_CURRENT_BINARY_DIR}/include")
  endif()
  if (HDF5_LIBRARY)
    if (NOT EXISTS ${HDF5_LIBRARY})
      message(FATAL_ERROR "Couldn't find HDF5 library at ${HDF5_LIBRARY}.")
    endif()
    message(STATUS "Using HDF5 library at ${HDF5_LIBRARY}.")
  else()
    set(HDF5_LIBRARY_DIR "${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")
    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
      set(HDF5_LIBRARY "${HDF5_LIBRARY_DIR}/libhdf5_debug.a")
    else()
      set(HDF5_LIBRARY "${HDF5_LIBRARY_DIR}/libhdf5.a")
    endif()
    message(STATUS "Building HDF5 library: ${HDF5_LIBRARY}.")
  endif()
  get_filename_component(HDF5_LIBRARY_DIR ${HDF5_LIBRARY} DIRECTORY)
  if (HDF5_HL_LIBRARY)
    if (NOT EXISTS ${HDF5_HL_LIBRARY})
      message(FATAL_ERROR "Couldn't find high-level HDF5 library at ${HDF5_HL_LIBRARY}.")
    endif()
    message(STATUS "Using high-level HDF5 library at ${HDF5_HL_LIBRARY}.")
  else()
    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
      set(HDF5_HL_LIBRARY "${HDF5_LIBRARY_DIR}/libhdf5_hl_debug.a")
    else()
      set(HDF5_HL_LIBRARY "${HDF5_LIBRARY_DIR}/libhdf5_hl.a")
    endif()
    message(STATUS "Building high-level HDF5 library: ${HDF5_HL_LIBRARY}.")
  endif()

  if (NETCDF_INCLUDE_DIR)
    if (NOT EXISTS ${NETCDF_INCLUDE_DIR})
      message(FATAL_ERROR "Couldn't find NetCDF include dir at ${NETCDF_INCLUDE_DIR}.")
    endif()
    message(STATUS "Using NetCDF include dir: ${NETCDF_INCLUDE_DIR}.")
  else()
    set(NETCDF_INCLUDE_DIR "${CMAKE_CURRENT_BINARY_DIR}/include")
  endif()
  if (NETCDF_LIBRARY)
    if (NOT EXISTS ${NETCDF_LIBRARY})
      message(FATAL_ERROR "Couldn't find NetCDF library at ${NETCDF_LIBRARY}.")
    endif()
    message(STATUS "Using NetCDF library at ${NETCDF_LIBRARY}.")
  else()
    set(NETCDF_LIBRARY_DIR "${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")
    set(NETCDF_LIBRARY "${NETCDF_LIBRARY_DIR}/libnetcdf.a")
    message(STATUS "Building NetCDF library: ${NETCDF_LIBRARY}.")
  endif()
  get_filename_component(NETCDF_LIBRARY_DIR ${NETCDF_LIBRARY} DIRECTORY)

  if (EKAT_INCLUDE_DIR)
    if (NOT EXISTS ${EKAT_INCLUDE_DIR})
      message(FATAL_ERROR "Couldn't find ekat include dir at ${EKAT_INCLUDE_DIR}.")
    endif()
    message(STATUS "Using ekat include dir: ${EKAT_INCLUDE_DIR}.")
  else()
    set(EKAT_INCLUDE_DIR "${CMAKE_CURRENT_BINARY_DIR}/include")
  endif()
  if (EKAT_LIBRARY)
    if (NOT EXISTS ${EKAT_LIBRARY})
      message(FATAL_ERROR "Couldn't find EKAT library at ${EKAT_LIBRARY}.")
    endif()
    message(STATUS "Using EKAT library at ${EKAT_LIBRARY}.")
  else()
    set(EKAT_LIBRARY_DIR "${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")
    set(EKAT_LIBRARY "${EKAT_LIBRARY_DIR}/libekat.a")
    message(STATUS "Building EKAT library: ${EKAT_LIBRARY}.")
  endif()
  get_filename_component(EKAT_LIBRARY_DIR ${EKAT_LIBRARY} DIRECTORY)

  if (YAMLCPP_INCLUDE_DIR)
    if (NOT EXISTS ${YAMLCPP_INCLUDE_DIR})
      message(FATAL_ERROR "Couldn't find yaml-cpp include dir at ${YAMLCPP_INCLUDE_DIR}.")
    endif()
    message(STATUS "Using yaml-cpp include dir: ${YAMLCPP_INCLUDE_DIR}.")
  else()
    set(YAMLCPP_INCLUDE_DIR "${CMAKE_CURRENT_BINARY_DIR}/include")
  endif()
  if (YAMLCPP_LIBRARY)
    if (NOT EXISTS ${YAMLCPP_LIBRARY})
      message(FATAL_ERROR "Couldn't find yaml-cpp library at ${YAMLCPP_LIBRARY}.")
    endif()
    message(STATUS "Using yaml-cpp library at ${YAMLCPP_LIBRARY}.")
  else()
    set(YAMLCPP_LIBRARY_DIR "${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")
    set(YAMLCPP_LIBRARY ${YAMLCPP_LIBRARY_DIR}/libyaml-cpp.a)
    message(STATUS "Building yaml-cpp library: ${YAMLCPP_LIBRARY}.")
  endif()
  get_filename_component(YAMLCPP_LIBRARY_DIR ${YAMLCPP_LIBRARY} DIRECTORY)

  if (OPENBLAS_INCLUDE_DIR)
    if (NOT EXISTS ${OPENBLAS_INCLUDE_DIR})
      message(FATAL_ERROR "Couldn't find OpenBLAS include dir at ${OPENBLAS_INCLUDE_DIR}.")
    endif()
    message(STATUS "Using OpenBLAS include dir: ${OPENBLAS_INCLUDE_DIR}.")
  else()
    set(OPENBLAS_INCLUDE_DIR "${CMAKE_CURRENT_BINARY_DIR}/include")
  endif()
  if (OPENBLAS_LIBRARY)
    if (NOT EXISTS ${OPENBLAS_LIBRARY})
      message(FATAL_ERROR "Couldn't find OpenBLAS library at ${OPENBLAS_LIBRARY}.")
    endif()
    message(STATUS "Using OpenBLAS library at ${OPENBLAS_LIBRARY}.")
  else()
    set(OPENBLAS_LIBRARY_DIR "${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")
    set(OPENBLAS_LIBRARY ${OPENBLAS_LIBRARY_DIR}/libopenblas.a)
    message(STATUS "Building OpenBLAS library: ${OPENBLAS_LIBRARY}.")
  endif()
  get_filename_component(OPENBLAS_LIBRARY_DIR ${OPENBLAS_LIBRARY} DIRECTORY)

  if (TINES_INCLUDE_DIR)
    if (NOT EXISTS ${TINES_INCLUDE_DIR})
      message(FATAL_ERROR "Couldn't find Tines include dir at ${TINES_INCLUDE_DIR}.")
    endif()
    message(STATUS "Using Tines include dir: ${TINES_INCLUDE_DIR}.")
  else()
    set(TINES_INCLUDE_DIR "${CMAKE_CURRENT_BINARY_DIR}/include")
  endif()
  if (TINES_LIBRARY)
    if (NOT EXISTS ${TINES_LIBRARY})
      message(FATAL_ERROR "Couldn't find Tines library at ${TINES_LIBRARY}.")
    endif()
    message(STATUS "Using Tines library at ${TINES_LIBRARY}.")
  else()
    set(TINES_LIBRARY_DIR "${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")
    set(TINES_LIBRARY ${TINES_LIBRARY_DIR}/libtines.a)
    message(STATUS "Building Tines library: ${TINES_LIBRARY}.")
  endif()
  get_filename_component(TINES_LIBRARY_DIR ${TINES_LIBRARY} DIRECTORY)

  if (TCHEM_INCLUDE_DIR)
    if (NOT EXISTS ${TCHEM_INCLUDE_DIR})
      message(FATAL_ERROR "Couldn't find TChem include dir at ${TCHEM_INCLUDE_DIR}.")
    endif()
    message(STATUS "Using TChem include dir: ${TCHEM_INCLUDE_DIR}.")
  else()
    set(TCHEM_INCLUDE_DIR "${CMAKE_CURRENT_BINARY_DIR}/include")
  endif()
  if (TCHEM_LIBRARY)
    if (NOT EXISTS ${TCHEM_LIBRARY})
      message(FATAL_ERROR "Couldn't find TChem library at ${TCHEM_LIBRARY}.")
    endif()
    message(STATUS "Using TChem library at ${TCHEM_LIBRARY}.")
  else()
    set(TCHEM_LIBRARY_DIR "${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")
    set(TCHEM_LIBRARY ${TCHEM_LIBRARY_DIR}/libtchem.a)
    message(STATUS "Building TChem library: ${TCHEM_LIBRARY}.")
  endif()
  get_filename_component(TCHEM_LIBRARY_DIR ${TCHEM_LIBRARY} DIRECTORY)

  if (APPLE)
    set(NEED_LAPACK FALSE)
  else()
    set(NEED_LAPACK TRUE)
  endif()

  # Certain tools (e.g. patch) require TMPDIR to be defined. If it is not,
  # we do so here.
  set(TMPDIR_VAR $ENV{TMPDIR})
  if (NOT TMPDIR_VAR)
    set(ENV{TMPDIR} "/tmp")
  endif()

endmacro()
