add_executable(haero_driver main.cpp driver.cpp)
target_link_libraries(haero_driver ${HAERO_LIBRARIES})
add_dependencies(haero_driver haero)
install(TARGETS haero_driver DESTINATION bin)

# Tests
#add_subdirectory(tests)
