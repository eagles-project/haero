add_executable(haero_driver main.cpp driver.cpp yaml_file.cpp)
target_link_libraries(haero_driver ${HAERO_DRIVER_LIBRARIES};${HAERO_LIBRARIES})
add_dependencies(haero_driver haero yaml_cpp)
install(TARGETS haero_driver DESTINATION bin)

# Tests
#add_subdirectory(tests)
