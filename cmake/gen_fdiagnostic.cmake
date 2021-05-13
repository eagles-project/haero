# This macro generates source files that implement a bridge between
# Haero's C++ interface and a Fortran module implementing a
# diagnostic function. Inputs:
# * CXX_CLASS_NAME, the name of a C++ class that will expose this process
# * PROCESS_ENUM_TYPE, a C++ enumerated type that defines the prognostic
#   process's role in the aerosol lifecycle
# * F90_MODULE_NAME, the name of a Fortran module that implements a prognostic
#   process
# Relevant source files are generated and/or appended to HAERO_DIAGNOSTICS, which
# is referenced in haero/diagnostics/CMakeLists.txt. This macro must be called
# only within that file.
function(faerosol_process CXX_CLASS_NAME PROCESS_ENUM_TYPE F90_MODULE_NAME)
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

function(faerosol_diagnostic CXX_CLASS_NAME PROCESS_ENUM_TYPE F90_MODULE_NAME)
 list(APPEND HAERO_DIAGNOSTICS ${CMAKE_CURRENT_SOURCE_DIR}/${F90_MODULE_NAME}.F90)
 set(CXX_HEADER ${CMAKE_CURRENT_BINARY_DIR}/${F90_MODULE_NAME}_fdiagnostic.hpp)
 configure_file(${CMAKE_CURRENT_SOURCE_DIR}/faerosol_diagnostic.hpp.in
                ${CXX_HEADER}
                @ONLY)

 set(CXX_SOURCE ${CMAKE_CURRENT_BINARY_DIR}/${F90_MODULE_NAME}_fdiagnostic.cpp)
 configure_file(${CMAKE_CURRENT_SOURCE_DIR}/faerosol_diagnostic.cpp.in
                ${CXX_SOURCE}
                @ONLY)
 list(APPEND HAERO_PROCESSES ${CXX_SOURCE})

 set(F90_BRIDGE ${CMAKE_CURRENT_BINARY_DIR}/${F90_MODULE_NAME}_bridge.F90)
 configure_file(${CMAKE_CURRENT_SOURCE_DIR}/diagnostic_bridge.F90.in
                ${F90_BRIDGE}
                @ONLY)
 list(APPEND HAERO_PROCESSES ${F90_BRIDGE})
endfunction()
