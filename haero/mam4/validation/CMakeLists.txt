# Here's a submodule that contains Skywalker inputs and Python scripts for
# processing outputs.
set(MAM_X_VALIDATION_DIR ${CMAKE_CURRENT_SOURCE_DIR}/mam_x_validation)

# This library has some common utilities used by Skywalker drivers in this
# directory tree.
add_library(validation validation.cpp)
add_dependencies(validation ext_libraries)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# These subdirectories contain Skywalker drivers for MAM4 parameterizations.
add_subdirectory(nucleation)
#add_subdirectory(gasaerexch)
