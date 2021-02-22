# Add source files for prognostic and diagnostic processes here.

# C++ process implementations go here. These are just names of C++ source
# files that implement aerosol prognostic and diagnostic processes.
set(CXX_PROCESS_SOURCES
#    mam4_nucleation_process.cpp
)

# Prepend processes/ to each source file and transmit HAERO_PROCESSES upward.
foreach(p_src ${CXX_PROCESS_SOURCES})
  list(APPEND HAERO_PROCESSES processes/${p_src})
endforeach()

if (HAERO_FORTRAN)
  include(gen_fprocess)
  include_directories(${CMAKE_CURRENT_BINARY_DIR})

  # List Fortran processes here. Declare a prognostic process with
  # the prognostic_fprocess macro, and a diagnostic process with the
  # diagnostic_fprocess macro. Both macros take the following arguments:
  # * The name of the C++ class to be generated for this process.
  # * The appropriate enumerated type that defines the role of the process
  #   within the aerosol life cycle (e.g. NucleationProcess for aerosol
  #   nucleation). For valid enum types, see ProcessType in haero/process.hpp.
  # * The name of the Fortran module that implements the process (without
  #   the .F90 suffix). The source code for this module must be located in a
  #   file named after the module, with an .F90 suffix. E.g.: mam4_nucleation
  #   refers to a Fortran module named mam4_nucleation whose source lives in
  #   mam4_nucleation.F90, within this directory.
  prognostic_fprocess(MAM4NucleationFProcess NucleationProcess mam4_nucleation)
  #diagnostic_fprocess(MAM4WaterUptakeFProcess WaterUptakeProcess mam4_water_uptake)

endif()

set(HAERO_PROCESSES ${HAERO_PROCESSES} PARENT_SCOPE)
