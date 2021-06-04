set(CHEM_DRIVER_INCDIRS ${CMAKE_CURRENT_BINARY_DIR})

include_directories(${HAERO_INCDIRS} ${CHEM_DRIVER_INCDIRS})

add_library(chem_standalone chemUtil.cpp
                            read_chem_input.cpp
                            chemDriver.cpp)

add_dependencies(chem_standalone haero yaml_cpp)

add_executable(chem_driver main.cpp)
target_link_libraries(chem_driver toy_problem;chem_standalone;yaml_cpp;${HAERO_LIBRARIES})

install(TARGETS chem_driver DESTINATION bin)
install(TARGETS chem_standalone DESTINATION lib)

# Tests
add_subdirectory(tests)
