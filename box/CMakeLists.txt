include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${PROJECT_BINARY_DIR}/haero)

add_executable(box main.cpp box_input.cpp box_input.F90 box_output.cpp)
target_link_libraries(box netcdf;${HAERO_LIBRARIES})

install(TARGETS box DESTINATION bin)

