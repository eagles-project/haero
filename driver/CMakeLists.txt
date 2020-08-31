add_executable(haero main.cpp driver.cpp)
target_link_libraries(haero ${HAERO_LIBRARIES})
install(TARGETS haero DESTINATION bin)

# Tests
#add_subdirectory(tests)
