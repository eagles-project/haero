include_directories(${CMAKE_CURRENT_BINARY_DIR} ${HAERO_INCDIRS} ${HAERO_DRIVER_INCDIRS})

add_executable(skywalker main.cpp parse_yaml.cpp)
target_link_libraries(skywalker yaml_cpp;${HAERO_LIBRARIES})

install(TARGETS skywalker DESTINATION bin)

