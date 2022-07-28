cmake_minimum_required (VERSION 3.12.0)

# Adjust CMake's module path.
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/ext/ekat/cmake/")

#---------
# Options
#---------
option(HAERO_ENABLE_CHEMISTRY "Enable TChem chemistry solver"              OFF)
option(HAERO_ENABLE_GPU       "Enable GPU support"                         OFF)
option(HAERO_ENABLE_MPI       "Enable MPI parallelism"                     ON)
option(HAERO_ENABLE_OPENMP    "Enable OpenMP threading (CPU only)"         OFF)
option(HAERO_ENABLE_MAM4BOX   "Enable MAM4 box model and related testing"  OFF)

# This option is only used in CI, where we have to use special sauce to get the
# submodules working with SSH. No mortal user should be concerned with this.
option(DISABLE_SUBMODULE_CHECKS "Skip checks for git submodules" OFF)
mark_as_advanced(FORCE DISABLE_SUBMODULE_CHECKS)

message(STATUS "Configuring with build type: ${CMAKE_BUILD_TYPE}")

# Check our pack size
if (NOT HAERO_PACK_SIZE)
  message(FATAL_ERROR "HAERO_PACK_SIZE is not specified! (Must be a positive integer).")
elseif(HAERO_PACK_SIZE LESS 1)
  message(FATAL_ERROR "Invalid HAERO_PACK_SIZE: ${HAERO_PACK_SIZE} (Must be a positive integer).")
endif()
message(STATUS "Configuring with pack size = ${HAERO_PACK_SIZE}")

# Configure the requested device (CPU or GPU).
if (NOT HAERO_DEVICE_ARCH)
  message(FATAL_ERROR "HAERO_DEVICE_ARCH was not provided!")
endif()

if(HAERO_ENABLE_GPU)
  message(STATUS "Building for GPU (${HAERO_DEVICE_ARCH})")
else() #CPU
  message(STATUS "Building for CPU (${HAERO_DEVICE_ARCH})")
endif()

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

enable_language(CXX)
enable_language(Fortran)

project (haero)

# Set all installation folders for third-party libraries, and figure out which
# ones have to be built for Haero.
include(set_up_platform)
set_up_platform()

# We use C++17 these days.
set(CMAKE_CXX_STANDARD 17)

message(STATUS "Generating project files in build directory: ${PROJECT_BINARY_DIR}")
message(STATUS "C++ compiler is ${CMAKE_CXX_COMPILER} (${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION})")
message(STATUS "C compiler is ${CMAKE_C_COMPILER} (${CMAKE_C_COMPILER_ID} ${CMAKE_C_COMPILER_VERSION})")
message(STATUS "Fortran compiler is ${CMAKE_Fortran_COMPILER} (${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION})")

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

# Fortran compiler flags
set(CMAKE_Fortran_FLAGS "-cpp")
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

# OpenMP support?
if (HAERO_ENABLE_OPENMP)
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
if (HAERO_ENABLE_GPU)
  set(CMAKE_CXX_COMPILER ${EKAT_NVCC_WRAPPER})
  # the esa_on_defaulted_function_ignored supresses a warning about KOKKOS markups of defaulted constructors/destructors
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --expt-extended-lambda --expt-relaxed-constexpr -Xcudafe --diag_suppress=esa_on_defaulted_function_ignored")
  # Avoid buggy experimental Kokkos CUDA stuff.
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DKOKKOS_CUDA_HALF_HPP_  -DKokkos_CXX_STANDARD:STRING=C++17")
else()
  # Make sure we're building with OpenMP where appropriate.
  if (NOT APPLE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
  endif()
endif()

# Record the libraries we've gathered so far as the "base" libraries that
# we need. We use this for efficiently linking unit tests.
set(HAERO_BASE_LIBRARIES ${HAERO_LIBRARIES})

# Include the binary directory in the header file search path,
# since it's where we place the third-party libraries.
include_directories("${PROJECT_BINARY_DIR}")
include_directories("${PROJECT_BINARY_DIR}/include")
include_directories(${HAERO_INCLUDE_DIRS})

include(GNUInstallDirs)
link_directories("${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")

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

# Testing
include(CTest)
enable_testing()

# Source directories.
include_directories("${PROJECT_SOURCE_DIR}")
add_subdirectory(haero)

if (HAERO_ENABLE_CHEMISTRY)
  message(STATUS "Building chem_driver.")
  include_directories("${TINES_INCLUDE_DIR}/tines")
  add_subdirectory(chem_driver)
endif()

# MAM4 box model (for testing ported parameterizations)
if (HAERO_ENABLE_MAM4BOX)
  message(STATUS "Enabling cross validation with MAM4 box model.")
  if (NOT MAM4BOX_DIR)
    message(FATAL_ERROR "MAM4 box model cross validation requested, but MAM4BOX_DIR is not set!")
  elseif(NOT EXISTS ${MAM4BOX_DIR})
    message(FATAL_ERROR "Invalid MAM4BOX_DIR: ${MAM4BOX_DIR}")
  endif()
  # Build the MAM4 box model's core library.
  set(MAM4BOX_CMAKE_OPTS -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}
                         -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                         -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                         -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                         -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
                         -DMAMBOX_PRECISION=${HAERO_PRECISION}
                         -DMAMBOX_ENABLE_DRIVER=OFF
                         -DMAMBOX_ENABLE_VALIDATION=OFF
                         -DMAMBOX_NUM_COLUMNS=1
                         -DMAMBOX_NUM_LEVELS=1
                         -DMAMBOX_NUM_SPECIES=35
                         -DMAMBOX_NUM_BLACK_CARBON_TAGS=1
                         -DMAMBOX_NUM_PRIMARY_ORGANIC_TAGS=1
                         -DMAMBOX_NUM_SECONDARY_ORGANIC_TAGS=1
                         -DMAMBOX_NUM_SECONDARY_ORGANIC_GAS_TAGS=1
                         -DMAMBOX_RAIN_EVAP_TO_COARSE_AERO=1)
  ExternalProject_Add(mam4box_proj
                      PREFIX ${CMAKE_CURRENT_BINARY_DIR}/mam4box
                      SOURCE_DIR ${MAM4BOX_DIR}
                      BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/mam4box
                      INSTALL_DIR ${PROJECT_BINARY_DIR}
                      CMAKE_ARGS ${MAM4BOX_CMAKE_OPTS}
                      LOG_CONFIGURE TRUE
                      BUILD_COMMAND ${MAKE} -j
                      LOG_BUILD TRUE
                      INSTALL_COMMAND ${MAKE} install
                      LOG_INSTALL TRUE)
  add_library(mam_core STATIC IMPORTED GLOBAL)
  set_target_properties(mam_core PROPERTIES IMPORTED_LOCATION ${PROJECT_BINARY_DIR}/lib/libmam_core.a)
  add_dependencies(mam_core mam4box_proj skywalker_f90)

  #----------------------------------------------
  # Build the MAM4 box model's Skywalker drivers
  #----------------------------------------------

  # Make sure the drivers can access the relevant Fortran modules and libraries.
  include_directories(${PROJECT_BINARY_DIR}/mam4box/core)
  #include_directories(${CMAKE_CURRENT_BINARY_DIR})
  set(MAMBOX_LIBRARIES skywalker_f90;skywalker;mam_core)

  # Copy the box model's standalone tests into its own source directory.
  file(COPY ${MAM4BOX_DIR}/standalone/tests/
       DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/mam4box_drivers)

  # Build the drivers.
  add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/mam4box_drivers)
endif()

# Formatting and format checking using clang-format.
find_program(CLANG_FORMAT clang-format)
if (NOT CLANG_FORMAT STREQUAL "CLANG_FORMAT-NOTFOUND")
  add_custom_target(format-cxx
    find ${PROJECT_SOURCE_DIR}/haero -name "*.[hc]pp" -exec ${CLANG_FORMAT} -i --style=google {} \+;
    VERBATIM
    COMMENT "Auto-formatting C++ code...")
  add_custom_target(format-cxx-check
    find ${PROJECT_SOURCE_DIR}/haero -name "*.[hc]pp" -exec ${CLANG_FORMAT} -n --Werror -ferror-limit=1 --style=google {} \+;
    VERBATIM
    COMMENT "Checking C++ formatting...")
endif()
