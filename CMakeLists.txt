# Minimum CMake version -- currently determined by HDF5.
cmake_minimum_required (VERSION 3.12.0)

# Adjust CMake's module path.
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/ext/ekat/cmake/")

set(ENABLE_FORTRAN ON)

# Check our pack size
if (NOT HAERO_PACK_SIZE)
  message(FATAL_ERROR "HAERO_PACK_SIZE is not specified! (Must be a positive integer).")
elseif(HAERO_PACK_SIZE LESS 1)
  message(FATAL_ERROR "Invalid HAERO_PACK_SIZE: ${HAERO_PACK_SIZE} (Must be a positive integer).")
elseif(HAERO_PACK_SIZE GREATER 1)
  set(ENABLE_FORTRAN OFF)
endif()
message(STATUS "Configuring with pack size = ${HAERO_PACK_SIZE}")

# This option is only used in CI, where we have to use special sauce to get the
# submodules working with SSH. No mortal user should be concerned with this.
option(HAERO_DISABLE_SUBMODULE_CHECKS "Skip checks for git submodules" OFF)
mark_as_advanced(FORCE HAERO_DISABLE_SUBMODULE_CHECKS)

# Check our device parameters
if (NOT HAERO_DEVICE)
  message(FATAL_ERROR "HAERO_DEVICE is not specified! (Must be CPU or CUDA).")
elseif(NOT HAERO_DEVICE STREQUAL CPU AND NOT HAERO_DEVICE STREQUAL CUDA)
  message(FATAL_ERROR "Invalid HAERO_DEVICE: ${HAERO_DEVICE} (Must be CPU or CUDA).")
elseif(HAERO_DEVICE STREQUAL CUDA)
  set(ENABLE_FORTRAN OFF)
  set(HAERO_USE_CUDA ON)
else() #CPU
  # If we're C++ only and not using a Mac, use OpenMP.
  if (NOT APPLE AND NOT ENABLE_FORTRAN)
    set(USE_OPENMP ON)
  else()
    set(USE_OPENMP OFF)
  endif()
endif()
message(STATUS "Configuring for Device: ${HAERO_DEVICE} (${HAERO_DEVICE_ARCH})")
message(STATUS "Configuring with build type: ${CMAKE_BUILD_TYPE}")

# Define default Kokkos View types
if (NOT HAERO_COLUMN_VIEW_TYPE)
  set(HAERO_COLUMN_VIEW_TYPE "DeviceType::view_1d<PackType>")
endif()
if (NOT HAERO_SPECIES_COLUMN_VIEW_TYPE)
  set(HAERO_SPECIES_COLUMN_VIEW_TYPE "DeviceType::view_2d<PackType>")
endif()
if (NOT HAERO_MODE_COLUMN_VIEW_TYPE)
  set(HAERO_MODE_COLUMN_VIEW_TYPE "DeviceType::view_2d<PackType>")
endif()

if (NOT ENABLE_FORTRAN)
  message(STATUS "================================================================")
  message(STATUS "  DISABLING FORTRAN PROCESSES: (GPU arch and/or pack size > 1)")
  message(STATUS "================================================================")
endif()

# We declare the project right away.
enable_language(CXX)
if (ENABLE_FORTRAN)
  set(HAERO_FORTRAN 1)
  enable_language(Fortran)
else()
  set(HAERO_FORTRAN 0)
endif()
project (haero)

# Set all installation folders for third-party libraries, and figure out which
# ones have to be built for Haero.
include(set_up_platform)
set_up_platform()

# We use C++14 these days.
set(CMAKE_CXX_STANDARD 14)

message(STATUS "Generating project files in build directory: ${PROJECT_BINARY_DIR}")
message(STATUS "C++ compiler is ${CMAKE_CXX_COMPILER} (${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION})")
message(STATUS "C compiler is ${CMAKE_C_COMPILER} (${CMAKE_C_COMPILER_ID} ${CMAKE_C_COMPILER_VERSION})")
if (HAERO_FORTRAN)
  message(STATUS "Fortran compiler is ${CMAKE_Fortran_COMPILER} (${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION})")
endif()

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
  set(CMAKE_Fortran_STANDARD_LIBRARIES "-lc++")
  message(STATUS "Configuring on Mac.")
elseif(LINUX)
  set(SYS_FLAGS "-DLINUX=1")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pthread")
  set(CMAKE_Fortran_STANDARD_LIBRARIES "-lstdc++ -no-pie")
  message(STATUS "Configuring on Linux.")
else()
  message(FATAL_ERROR "Unknown system type! Currently only Mac and Linux are supported")
endif ()

# C++ compiler flags.
set(CMAKE_CXX_FLAGS_RELEASE "-O2 -DNDEBUG")
set(CMAKE_CXX_FLAGS_RelWithDebInfo "-g -O2 -DNDEBUG")
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

if (HAERO_FORTRAN)
  # Fortran compiler flags
  set(CMAKE_Fortran_FLAGS "-cpp -Dc_real=c_${HAERO_REAL_TYPE}")
  if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fbacktrace")
    endif()
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wall -ffree-form -ffree-line-length-none")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wno-unused-dummy-argument")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Werror=tabs")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wno-maybe-uninitialized")
    if (APPLE)
      # Mac-specific stuff goes here.
    endif()
  elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -traceback")
    endif()
  endif()
endif()

# Report the installation prefix.
message(STATUS "Installation prefix is ${CMAKE_INSTALL_PREFIX}")

# Basic libraries to be linked in.
set(HAERO_LIBRARIES m)

# Process any additional linker flags supplied with EXTRA_LDFLAGS.
if (HAERO_EXTRA_LDFLAGS)
  include(parse_extra_ldflags)
  parse_extra_ldflags(${HAERO_EXTRA_LDFLAGS})
endif()

# Figure out MPI.
if (HAVE_MPI)
  message(STATUS "MPI is enabled")
  set(HAVE_MPI 1) # Normalize.

  # NOTE: Disable C++ bindings for MPI, since they have never worked for anyone.
  set(NO_MPI_CXX_FLAGS "-DMPICH_SKIP_MPICXX -UHAVE_MPI_CPP -DLAM_WANT_MPI2CPP=0 -DLAM_BUILDING=1 -DOMPI_WANT_CXX_BINDINGS=0 -DOMPI_BUILDING=1")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${NO_MPI_CXX_FLAGS}")
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
  set(HAVE_MPI 0) # Normalize.
endif()

# OpenMP support?
if (USE_OPENMP)
  find_package(OpenMP QUIET)
  if (OPENMP_FOUND)
    message(STATUS "OpenMP is enabled")
    set(HAVE_OPENMP 1)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  else()
    message(STATUS "OpenMP is not supported")
    set(HAVE_OPENMP 0)
  endif()
else()
  message(STATUS "OpenMP is disabled")
  set(HAVE_OPENMP 0)
endif()

# Other third-party libraries.
add_subdirectory(ext)

# If we're building on a GPU, set the C++ compiler appropriately.
if (HAERO_DEVICE STREQUAL CUDA)
  set(CMAKE_CXX_COMPILER "${EKAT_LIBRARY_DIR}/../bin/ekat_mpicxx")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --expt-extended-lambda --expt-relaxed-constexpr")
  # Avoid buggy experimental Kokkos CUDA stuff.
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DKOKKOS_CUDA_HALF_HPP_")

endif()

# Record the libraries we've gathered so far as the "base" libraries that
# our box model needs. We use this for efficiently linking unit tests.
set(HAERO_BASE_LIBRARIES ${HAERO_LIBRARIES})

# Include the binary directory in the header file search path,
# since it's where we place the third-party libraries.
include_directories("${PROJECT_BINARY_DIR}")
include_directories("${PROJECT_BINARY_DIR}/include")
include_directories("${TINES_INCLUDE_DIR}/tines")

include(GNUInstallDirs)
link_directories("${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")
include_directories(${HAERO_INCDIRS})

# Support for valgrind -- Linux only.
if (LINUX)
  find_package(Valgrind QUIET)
  if (VALGRIND_FOUND)
    set(VALGRIND_FOUND 1) # regularize this value
    include_directories(${VALGRIND_INCLUDE_DIR})
    set(MEMORYCHECK_COMMAND ${VALGRIND_PROGRAM})
    # Add "--gen-suppressions=all" to MEMORYCHECK_COMMAND_OPTIONS to generate
    # suppressions for Valgrind's false positives. The suppressions show up
    # right in the MemoryChecker.*.log files.
    set(MEMORYCHECK_COMMAND_OPTIONS "--leak-check=full --show-leak-kinds=all --errors-for-leak-kinds=definite,possible --track-origins=yes --error-exitcode=1 --trace-children=yes --suppressions=${PROJECT_SOURCE_DIR}/tools/valgrind/scasm.supp" CACHE STRING "Options passed to Valgrind." FORCE)

    # make memcheck target
    add_custom_target(memcheck ctest -T memcheck -j ${NUMBER_OF_CORES} USES_TERMINAL)
  else()
    set(VALGRIND_FOUND 0)
  endif()
else()
  # Valgrind doesn't work on Macs.
  set(VALGRIND_FOUND 0)
endif()

# Unit testing.
include(CTest)
enable_testing()

# Source directories.
include_directories("${PROJECT_SOURCE_DIR}")
add_subdirectory(haero)

if (HAERO_ENABLE_DRIVER)
  message(STATUS "Building standalone haero_driver.")
  add_subdirectory(driver)
endif()

#if (HAERO_ENABLE_CHEM_DRIVER)
#  message(STATUS "Building chem_driver.")
#  add_subdirectory(chem_driver)
#endif()

# Validation tools (currently offered only for CPU builds with pack size 1,
# because we're using it with Fortran.)
if(NOT HAERO_USE_CUDA AND HAERO_PACK_SIZE EQUAL 1)
  # Aerosol cross-validation materials
  message(STATUS "Building validation programs...")
  add_subdirectory(validation)
endif()

# Formatting and format checking using clang-format.
find_program(CLANG_FORMAT clang-format)
if (NOT CLANG_FORMAT STREQUAL "CLANG_FORMAT-NOTFOUND")
  add_custom_target(format
    find ${PROJECT_SOURCE_DIR}/haero -name "*.[hc]pp" -exec ${CLANG_FORMAT} -i --style=google {} \+;
    COMMAND find ${PROJECT_SOURCE_DIR}/driver -name "*.[hc]pp" -exec ${CLANG_FORMAT} -i --style=google {} \+;
    COMMAND find ${PROJECT_SOURCE_DIR}/validation -name "*.[hc]pp" -exec ${CLANG_FORMAT} -i --style=google {} \+;
    VERBATIM
    COMMENT "Auto-formatting C++ code...")
  add_custom_target(format-check
    find ${PROJECT_SOURCE_DIR}/haero -name "*.[hc]pp" -exec ${CLANG_FORMAT} -n --Werror -ferror-limit=1 --style=google {} \+;
    COMMAND find ${PROJECT_SOURCE_DIR}/driver -name "*.[hc]pp" -exec ${CLANG_FORMAT} -n --Werror -ferror-limit=1 --style=google {} \+;
    COMMAND find ${PROJECT_SOURCE_DIR}/validation -name "*.[hc]pp" -exec ${CLANG_FORMAT} -n --Werror -ferror-limit=1 --style=google {} \+;
    VERBATIM
    COMMENT "Checking C++ formatting...")
endif()
