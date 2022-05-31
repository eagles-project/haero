set(MAM_X_VALIDATION_DIR ${CMAKE_CURRENT_SOURCE_DIR}/mam_x_validation)

# Here are C++ and Fortran libraries with common utilities used by validation
# drivers in the ѕubdirectories below.
add_library(validation validation.cpp)
add_dependencies(validation ext_libraries)
add_library(validation_f90 validation.F90)
add_dependencies(validation_f90 ext_libraries)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(nucleation)
add_subdirectory(gasaerexch)
