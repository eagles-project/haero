# This macro generates source files that implement a bridge between
# Haero's C++ interface and a Fortran module implementing a prognostic
# process. Inputs:
# * CXX_CLASS_NAME, the name of a C++ class that will expose this process
# * F90_MODULE_NAME, the name of a Fortran module that implements a prognostic
#   process
# Relevant source files are generated and/or appended to HAERO_PROCESSES, which
# is referenced in haero/processes/CMakeLists.txt. This macro must be called
# only within that file.
function(faerosol_process CXX_CLASS_NAME F90_MODULE_NAME)
  list(APPEND HAERO_PROCESSES ${CMAKE_CURRENT_SOURCE_DIR}/${F90_MODULE_NAME}.F90)
  set(CXX_HEADER ${CMAKE_CURRENT_BINARY_DIR}/${F90_MODULE_NAME}_fprocess.hpp)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/faerosol_process.hpp.in
                 ${CXX_HEADER}
                 @ONLY)

  set(CXX_SOURCE ${CMAKE_CURRENT_BINARY_DIR}/${F90_MODULE_NAME}_fprocess.cpp)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/faerosol_process.cpp.in
                 ${CXX_SOURCE}
                 @ONLY)
  list(APPEND HAERO_PROCESSES ${CXX_SOURCE})

  set(F90_BRIDGE ${CMAKE_CURRENT_BINARY_DIR}/${F90_MODULE_NAME}_bridge.F90)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/faerosol_process_bridge.F90.in
                 ${F90_BRIDGE}
                 @ONLY)
  list(APPEND HAERO_PROCESSES ${F90_BRIDGE})
  set(HAERO_PROCESSES ${HAERO_PROCESSES} PARENT_SCOPE)
endfunction()

