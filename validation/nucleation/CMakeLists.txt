add_executable(validate_nucleation validate.cpp)
target_link_libraries(validate_nucleation skywalker;yaml_cpp;${HAERO_LIBRARIES})
