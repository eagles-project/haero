include_directories(${CMAKE_CURRENT_BINARY_DIR} ${HAERO_INCDIRS} ${HAERO_DRIVER_INCDIRS})

add_library(haero_standalone host_params.cpp
                             host_state.cpp
#                             dyn_column.cpp
                             read_yaml.cpp
                             driver.cpp
                             ncwriter.cpp
                             ncreader.cpp)
add_dependencies(haero_standalone haero yaml_cpp netcdf)

add_executable(haero_driver main.cpp)
target_link_libraries(haero_driver haero_standalone;${HAERO_DRIVER_LIBRARIES};${HAERO_LIBRARIES})

install(TARGETS haero_driver DESTINATION bin)
install(TARGETS haero_standalone DESTINATION lib)

# Tests
add_subdirectory(tests)
