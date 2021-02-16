# Add source files for prognostic and diagnostic processes here.

# C++ process implementations go here.
#set(CXX_PROCESS_SOURCES
#    mam4_nucleation_process.cpp
#)

# Fortran process implementations (including related C++ files) go here.
set(FORTRAN_PROCESS_SOURCES
    mam4_nucleation_fprocess.cpp
    mam4_nucleation_fprocess.F90
)

# Prepend processes/ to each source file and transmit HAERO_PROCESSES upward.
foreach(p_src ${CXX_PROCESS_SOURCES})
  list(APPEND HAERO_PROCESSES processes/${p_src})
endforeach()

if (HAERO_FORTRAN)
  foreach(p_src ${FORTRAN_PROCESS_SOURCES})
    list(APPEND HAERO_PROCESSES processes/${p_src})
  endforeach()
endif()

set(HAERO_PROCESSES ${HAERO_PROCESSES} PARENT_SCOPE)
