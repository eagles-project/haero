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

  # Configure EKAT (whether built externally or internally).
  if (EKAT_SOURCE_DIR)
    if (NOT EXISTS ${EKAT_SOURCE_DIR})
      message(FATAL_ERROR "Invalid EKAT source dir: ${EKAT_SOURCE_DIR}.")
    elseif()
      if (NOT EXISTS ${EKAT_SOURCE_DIR}/src)
        message(FATAL_ERROR "EKAT source dir (${EKAT_SOURCE_DIR}) has no src/ ѕubdirectory!")
      endif()
    endif()
    if (NOT EKAT_BINARY_DIR)
      message(FATAL_ERROR "EKAT source dir was given, but binary dir was not!")
    elseif (NOT EXISTS ${EKAT_BINARY_DIR})
      message(FATAL_ERROR "Invalid EKAT binary dir: ${EKAT_BINARY_DIR}.")
    elseif (NOT EXISTS ${EKAT_BINARY_DIR}/src/ekat/libekat.a)
      message(FATAL_ERROR "EKAT binary dir (${EKAT_BINARY_DIR}) does not contain an EKAT library!")
    endif()
    message(STATUS "Using pre-built EKAT library in ${EKAT_BINARY_DIR}.")
    set(HAERO_BUILDS_EKAT OFF)
  else()
    set(EKAT_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/ext/ekat")
    set(EKAT_BINARY_DIR "${CMAKE_CURRENT_BINARY_DIR}/externals/ekat")
    message(STATUS "Building EKAT internally.")
    set(HAERO_BUILDS_EKAT ON)
  endif()

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
