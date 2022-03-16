# Add source files for prognostic processes here.

# C++ process implementations go here. These are just names of C++ source
# files that implement aerosol prognostic and diagnostic processes.
set(CXX_PROCESS_SOURCES
    mam_nucleation_process.cpp
    mam_gasaerexch_process.cpp
    mam_calcsize_process.cpp
    mam_rename_process.cpp
    simple_nucleation_process.cpp
    chemistry_process.cpp
)

# This is a transitional process used for debugging on the host
if(HAERO_DEVICE STREQUAL CPU)
  list(APPEND CXX_PROCESS_SOURCES mam_calcsize_hostcxx_process.cpp)
endif()

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
  # * The name of the Fortran module that implements the process (without
  #   the .F90 suffix). The source code for this module must be located in a
  #   file named after the module, with an .F90 suffix. E.g.: mam_nucleation
  #   refers to a Fortran module named mam_nucleation whose source lives in
  #   mam_nucleation.F90, within this directory.
  faerosol_process(MAMGasAerosolExchangeFProcess mam_gasaerexch)
  faerosol_process(MAMNucleationFProcess mam_nucleation)
  # TODO: Temporarily disabled for Fortran modal_aero_config refactor -JNJ
  faerosol_process(MAMCalcsizeFProcess mam_calcsize)
  faerosol_process(MAMRenameFProcess mam_rename)

  # Supporting modules get appended to HAERO_FORTRAN_SOURCES here.
  foreach(f90_src vehkamaki2002.F90 merikanto2007.F90 wang2008.F90
          kerminen2002.F90)
    list(APPEND HAERO_FORTRAN_SOURCES processes/${f90_src})
  endforeach()
  set(HAERO_FORTRAN_SOURCES ${HAERO_FORTRAN_SOURCES} PARENT_SCOPE)

endif()

set(HAERO_PROCESSES ${HAERO_PROCESSES} PARENT_SCOPE)
