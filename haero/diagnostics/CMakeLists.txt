# Add source files for diagnostic functions here.

# C++ diagnostic implementations go here.  These are the names of the
# c++ source files that implement diagnostic functions
set (CXX_DIAGNOSTIC_SOURCES
  kohler_solve_diagnostic.cpp
  mode_dry_particle_volume_diagnostic.cpp
  modal_mean_hygroscopicity_diagnostic.cpp
  modal_wet_radius_diagnostic.cpp
)

# prepend diagnostics/ to each source file and transmit HAERO_DIAGNOSTICS upward
foreach(d_src ${CXX_DIAGNOSTIC_SOURCES})
  list(APPEND HAERO_DIAGNOSTICS diagnostics/${d_src})
endforeach()

if (HAERO_FORTRAN)
  include_directories(${CMAKE_CURRENT_BINARY_DIR})
  # TODO: declare fortran diagnostics here.
endif()

set(HAERO_DIAGNOSTICS ${HAERO_DIAGNOSTICS} PARENT_SCOPE)
