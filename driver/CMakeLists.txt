include_directories(${CMAKE_CURRENT_BINARY_DIR})

add_library(haero_standalone atmosphere.cpp
                             phys_column.cpp
                             read_yaml.cpp
                             driver.cpp)
add_dependencies(haero_standalone haero)

add_executable(haero_driver main.cpp)
                            dyn_column.cpp)

add_executable(haero_driver main.cpp driver.cpp read_yaml.cpp)
target_link_libraries(haero_driver haero_standalone;${HAERO_DRIVER_LIBRARIES};${HAERO_LIBRARIES})

install(TARGETS haero_driver DESTINATION bin)
install(TARGETS haero_standalone DESTINATION lib)

# Tests
add_subdirectory(tests)
