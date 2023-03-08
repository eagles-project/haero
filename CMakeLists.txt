cmake_minimum_required (VERSION 3.12.0)

# Adjust CMake's module path.
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/ext/ekat/cmake/")

#---------
# Options
#---------
option(HAERO_ENABLE_GPU       "Enable GPU support"                         OFF)
option(HAERO_ENABLE_MPI       "Enable MPI parallelism"                     ON)

# This option is only used in CI, where we have to use special sauce to get the
# submodules working with SSH. No mortal user should be concerned with this.
option(DISABLE_SUBMODULE_CHECKS "Skip checks for git submodules" OFF)
mark_as_advanced(FORCE DISABLE_SUBMODULE_CHECKS)

message(STATUS "Configuring with build type: ${CMAKE_BUILD_TYPE}")

if(HAERO_ENABLE_GPU)
  message(STATUS "Building for GPU")
else() #CPU
  message(STATUS "Building for CPU")
endif()

# Blessed version of clang-format.
set(CLANG_FORMAT_VERSION 14)

enable_language(C)   # needed so ekat can detect MPI
enable_language(CXX)

project (haero)

# Set all installation folders for third-party libraries, and figure out which
# ones have to be built for Haero.
include(set_up_platform)
set_up_platform()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_C_STANDARD 11)

message(STATUS "Generating project files in build directory: ${PROJECT_BINARY_DIR}")
message(STATUS "C++ compiler is ${CMAKE_CXX_COMPILER} (${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION})")
message(STATUS "C compiler is ${CMAKE_C_COMPILER} (${CMAKE_C_COMPILER_ID} ${CMAKE_C_COMPILER_VERSION})")

# Version numbers.
set (HAERO_MAJOR_VERSION 0)
set (HAERO_MINOR_VERSION 1)
set (HAERO_PATCH_VERSION 0)

set (HAERO_VERSION "${HAERO_MAJOR_VERSION}.${HAERO_MINOR_VERSION}.${HAERO_PATCH_VERSION}")
message(STATUS "Configuring haero v${HAERO_VERSION}")

# Precision of floating point numbers.
if (HAERO_PRECISION STREQUAL "single")
  set(HAERO_REAL_TYPE "float")
  set(HAERO_REAL_KIND "sp")
  set(HAERO_DOUBLE_PRECISION 0)
elseif(HAERO_PRECISION STREQUAL "double")
  set(HAERO_REAL_TYPE "double")
  set(HAERO_REAL_KIND "dp")
  set(HAERO_DOUBLE_PRECISION 1)
elseif(HAERO_PRECISION)
  message(FATAL_ERROR "Invalid HAERO_PRECISION: ${HAERO_PRECISION} (use 'single' or 'double')")
else()
  message(FATAL_ERROR "HAERO_PRECISION not set (use 'single' or 'double')")
endif()
message(STATUS "Using ${HAERO_PRECISION} precision floating point numbers")

# We build static libraries only.
set(BUILD_SHARED_LIBS OFF)

# Figure out the system type.
if (APPLE)
  set(SYS_FLAGS "-DAPPLE=1")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -framework Accelerate")
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -framework Accelerate")
  message(STATUS "Configuring on Mac.")
elseif(LINUX)
  set(SYS_FLAGS "-DLINUX=1")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pthread")
  message(STATUS "Configuring on Linux.")
else()
  message(FATAL_ERROR "Unknown system type! Currently only Mac and Linux are supported")
endif ()

# C++ compiler flags.
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-sign-compare")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fmax-errors=10")
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas")
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -mkl -static-intel")
  message(STATUS "Using Intel compilers and MKL linear algebra libraries")
endif()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SYS_FLAGS}")

# C compiler flags
set(CMAKE_C_FLAGS "-fPIC")
if (CMAKE_C_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-format-truncation -Wfloat-equal")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-unknown-pragmas")
  if (HAVE_DOUBLE_PRECISION)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wfloat-conversion")
  endif()

  if (LINUX)
    # Pass some more needed flags to the compiler.
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pthread")
  endif()
elseif (CMAKE_C_COMPILER_ID MATCHES "Clang")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-unknown-pragmas")

  # Clang proper behaves differently from Apple Clang, and requires
  # more suppressions.
  if (NOT CMAKE_C_COMPILER_ID MATCHES "AppleClang")
    # tolower() is a recursive macro, which causes issues.
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-disabled-macro-expansion")
  endif()
elseif (CMAKE_C_COMPILER_ID MATCHES "Intel")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall")
endif()
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${SYS_FLAGS}")

# Report the installation prefix.
message(STATUS "Installation prefix is ${CMAKE_INSTALL_PREFIX}")

# Basic libraries to be linked in.
set(HAERO_LIBRARIES m)

# Figure out MPI.
if (HAERO_ENABLE_MPI)
  message(STATUS "MPI is enabled")

  find_package(MPI REQUIRED COMPONENTS C)
  include(EkatMpiUtils)
  DisableMpiCxxBindings()

  include_directories(${MPI_C_INCLUDE_DIRS})
  set(HAERO_LIBRARIES ${MPI_C_LIBRARIES};${HAERO_LIBRARIES})

  if (HAERO_MPI_EXEC)
    message(STATUS "MPI exec: ${HAERO_MPI_EXEC}")
    if (NOT HAERO_MPI_NP_FLAG)
      message(FATAL_ERROR "HAERO_MPI_EXEC was given, but HAERO_MPI_NP_FLAG was not!")
    endif()
    message(STATUS "MPI np flag: ${HAERO_MPI_NP_FLAG}")
    if (HAERO_MPI_EXTRA_FLAGS)
      message(STATUS "MPI extra flags: ${HAERO_MPI_EXTRA_FLAGS}")
    endif()
  endif()
else()
  message(STATUS "MPI is disabled")
endif()

# Enable OpenMP for CPU backends if available.
if (NOT HAERO_ENABLE_GPU AND NOT APPLE)
  find_package(OpenMP QUIET)
  if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  endif()
endif()

# Other third-party libraries.
add_subdirectory(ext)

# If we're building on a GPU, set the C++ compiler appropriately.
if (HAERO_ENABLE_GPU)
  set(CMAKE_CXX_COMPILER ${EKAT_NVCC_WRAPPER})
  # the esa_on_defaulted_function_ignored supresses a warning about KOKKOS markups of defaulted constructors/destructors
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --expt-extended-lambda --expt-relaxed-constexpr -Xcudafe --diag_suppress=esa_on_defaulted_function_ignored")
  # Avoid buggy experimental Kokkos CUDA stuff.
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DKOKKOS_CUDA_HALF_HPP_  -DKokkos_CXX_STANDARD=C++17")
endif()

# Record the libraries we've gathered so far as the "base" libraries that
# we need. We use this for efficiently linking unit tests.
set(HAERO_BASE_LIBRARIES ${HAERO_LIBRARIES})

# Include the binary directory in the header file search path,
# since it's where we place the third-party libraries.
include_directories("${PROJECT_BINARY_DIR}")
include_directories("${PROJECT_BINARY_DIR}/include")
include_directories(${HAERO_INCLUDE_DIRS})
# Add external include files as SYSTEM to turn off warnings from these
include_directories(SYSTEM ${HAERO_EXT_INCLUDE_DIRS})

include(GNUInstallDirs)
link_directories("${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")

# Testing
include(CTest)
enable_testing()

# Source directories.
include_directories("${PROJECT_SOURCE_DIR}")
add_subdirectory(haero)

# Generate haero.cmake (for downstream aerosol packages) and install it.
set(HAERO_TPL_IMPORTED_LOCATIONS "") # <-- generated segment of haero.cmake
foreach (lib kokkoscore kokkoscontainers kokkossimd spdlog yaml-cpp ekat ekat_test_main ekat_test_session haero)
  string(APPEND HAERO_TPL_IMPORTED_LOCATIONS "add_library(${lib} STATIC IMPORTED GLOBAL)\n")
  set(suffix "")
  if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    # It seems to be a trend these days to append 'd' to your library in
    # Debug builds. :-/
    if (${lib} MATCHES "spdlog" OR ${lib} MATCHES "yaml-cpp")
      set(suffix "d")
    endif()
  endif()
  get_target_property(imported_loc ${lib} IMPORTED_LOCATION)
  if (imported_loc MATCHES "-NOTFOUND")
    set(imported_loc ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/lib${lib}${suffix}.a)
  endif()
  string(APPEND HAERO_TPL_IMPORTED_LOCATIONS "set_target_properties(${lib} PROPERTIES IMPORTED_LOCATION ${imported_loc})\n")
endforeach()
configure_file(
  ${PROJECT_SOURCE_DIR}/cmake/haero.cmake.in
  ${PROJECT_BINARY_DIR}/share/haero.cmake
  @ONLY
)
install(DIRECTORY ${PROJECT_BINARY_DIR}/share DESTINATION .)

# Support for valgrind -- Linux only.
if (LINUX)
  find_program(VALGRIND_EXE valgrind)
  if (NOT VALGRIND_EXE MATCHES "-NOTFOUND")
    set(VALGRIND_FOUND 1) # regularize this value
    set(MEMORYCHECK_COMMAND ${VALGRIND_EXE})
    # Add "--gen-suppressions=all" to MEMORYCHECK_COMMAND_OPTIONS to generate
    # suppressions for Valgrind's false positives. The suppressions show up
    # right in the MemoryChecker.*.log files.
    set(MEMORYCHECK_COMMAND_OPTIONS "--leak-check=full --show-leak-kinds=all --errors-for-leak-kinds=definite,possible --track-origins=yes --error-exitcode=1 --trace-children=yes" CACHE STRING "Options passed to Valgrind." FORCE)

    # make memcheck target
    add_custom_target(memcheck ctest -T memcheck -j ${NUMBER_OF_CORES} USES_TERMINAL)
  else()
    set(VALGRIND_FOUND 0)
  endif()
else()
  # Valgrind doesn't work on Macs.
  set(VALGRIND_FOUND 0)
endif()

# Formatting and format checking using clang-format.
find_program(CLANG_FORMAT clang-format)
if (NOT CLANG_FORMAT STREQUAL "CLANG_FORMAT-NOTFOUND")
  # Is this the blessed version? If not, we create targets that warn the user
  # to obtain the right version.
  execute_process(COMMAND clang-format --version
    OUTPUT_VARIABLE CF_VERSION)
  string(STRIP ${CF_VERSION} CF_VERSION)
  if (NOT ${CF_VERSION} MATCHES ${CLANG_FORMAT_VERSION})
    add_custom_target(format-cxx
      echo "You have clang-format version ${CF_VERSION}, but ${CLANG_FORMAT_VERSION} is required."
      "Please make sure this version appears in your path and rerun config.sh.")
    add_custom_target(format-cxx-check
      echo "You have clang-format version ${CF_VERSION}, but ${CLANG_FORMAT_VERSION} is required."
      "Please make sure this version appears in your path and rerun config.sh.")
  else()
    add_custom_target(format-cxx
      find ${PROJECT_SOURCE_DIR}/haero -name "*.[hc]pp" -exec ${CLANG_FORMAT} -i {} \+;
      VERBATIM
      COMMENT "Auto-formatting C++ code...")
    add_custom_target(format-cxx-check
      find ${PROJECT_SOURCE_DIR}/haero -name "*.[hc]pp" -exec ${CLANG_FORMAT} -n --Werror -ferror-limit=1 {} \+;
      VERBATIM
      COMMENT "Checking C++ formatting...")
  endif()

endif()
