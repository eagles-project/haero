add_executable(haero_driver main.cpp driver.cpp)
target_link_libraries(haero_driver ${HAERO_LIBRARIES})
install(TARGETS haero_driver DESTINATION bin)

# Tests
#add_subdirectory(tests)
