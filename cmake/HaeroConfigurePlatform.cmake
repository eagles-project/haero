# This macro identifies compilers and third-party library needs
# for particular hosts.
macro(HaeroConfigurePlatform)

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

  # Haero can be configured either as a standalone library or as part of a
  # larger project.
  #
  # If built as a standalone library, it can either use an
  # existing installation of ekat or build its own. This behavior is governed
  # by the HAERO_BUILDS_EKAT flag (if ON, Haero builds its own ekat, if OFF,
  # it uses an existing one).
  set(HAERO_BUILDS_EKAT OFF)

  # If built as part of a larger project, an ekat target must exist within
  # the CMake build system already. Haero thus uses the existence of an ekat
  # target to determine whether it's part of a larger project
  set(HAERO_STANDALONE ON)

  if (TARGET ekat) # we're part of a larger project
    message(STATUS "Building Haero within another project.")
    set(HAERO_STANDALONE OFF)
  else() # we're in standalone mode
    message(STATUS "Building Haero in standalone mode.")
    if (EKAT_SOURCE_DIR) # ekat paths specified
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
    else() # no specified paths to ekat
      find_package(ekat) # can we find it on the system?
      if (NOT ekat_FOUND) # no... build it ourselves
        set(EKAT_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/ext/ekat")
        set(EKAT_BINARY_DIR "${CMAKE_CURRENT_BINARY_DIR}/externals/ekat")
        message(STATUS "Building EKAT internally.")
        set(HAERO_BUILDS_EKAT ON)
      endif()
    endif()
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
