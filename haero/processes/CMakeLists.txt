# Add source files for prognostic and diagnostic processes here.
set(PROCESS_SOURCES
#    mam4_water_uptake.cpp
#    mam4_water_uptake.f90
)

# Prepend processes/ to each source file and transmit HAERO_PROCESSES upward.
foreach(p_src ${PROCESS_SOURCES})
  list(APPEND HAERO_PROCESSES processes/${p_src})
endforeach()
set(HAERO_PROCESSES ${HAERO_PROCESSES} PARENT_SCOPE)
