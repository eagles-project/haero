# Add source files for prognostic processes here.

# C++ process implementations go here. These are just names of C++ source
# files that implement aerosol prognostic and diagnostic processes.
set(CXX_PROCESS_SOURCES
    mam_nucleation_process.cpp
    mam_calcsize_process.cpp
    mam_rename_subarea_process.cpp
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
  #   file named after the module, with an .F90 suffix. E.g.: mam_nucleation
  #   refers to a Fortran module named mam_nucleation whose source lives in
  #   mam_nucleation.F90, within this directory.
  faerosol_process(MAMNucleationFProcess NucleationProcess mam_nucleation)
  faerosol_process(MAMCalcsizeFProcess CalcsizeProcess mam_calcsize)
  faerosol_process(MAMResizeSubareaFProcess RenameSubareaProcess mam_rename_subarea)

endif()

set(HAERO_PROCESSES ${HAERO_PROCESSES} PARENT_SCOPE)
