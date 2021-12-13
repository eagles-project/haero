set(MAM_X_VALIDATION_DIR ${CMAKE_CURRENT_SOURCE_DIR}/mam_x_validation)

# Here's a library of common utilities used by validation drivers in the
# ѕubdirectories below.
add_library(validation validation.cpp)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(nucleation)
