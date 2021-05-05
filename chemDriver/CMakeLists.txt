include_directories(${CMAKE_CURRENT_BINARY_DIR} ${HAERO_INCDIRS} ${CHEM_DRIVER_INCDIRS})

# add_library(chem_standalone host_params.cpp
#                             host_dynamics.cpp
#                             host_dynamics_tends.cpp
#                             read_yaml.cpp
#                             driver.cpp
#                             ncwriter.cpp
#                             ncreader.cpp)
add_library(chem_standalone chemUtil.cpp
                            read_chem_input.cpp
                            chemDriver.cpp)
# add_dependencies(chem_standalone haero yaml_cpp)
add_dependencies(chem_standalone haero yaml_cpp)

add_executable(chem_driver main.cpp)
# target_link_libraries(chem_driver chem_standalone;${CHEM_DRIVER_LIBRARIES};${HAERO_LIBRARIES})
target_link_libraries(chem_driver chem_standalone;yaml_cpp;${HAERO_LIBRARIES})

install(TARGETS chem_driver DESTINATION bin)
install(TARGETS chem_standalone DESTINATION lib)

# Tests
add_subdirectory(tests)
