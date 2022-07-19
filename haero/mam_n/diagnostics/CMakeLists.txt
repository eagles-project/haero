# Add source files for diagnostic functions here.

# C++ source files for diagnostic functions go here. Many diagnostic functions
# can be defined in headers, which don't need to appear here.
set (CXX_DIAGNOSTIC_SOURCES
  kohler_solve.cpp
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
