set(CHEM_DRIVER_INCDIRS ${CMAKE_CURRENT_BINARY_DIR})

include_directories(${HAERO_INCDIRS} ${CHEM_DRIVER_INCDIRS})

add_library(chem_standalone chem_driver.cpp)

add_dependencies(chem_standalone haero yaml_cpp)

install(TARGETS chem_standalone DESTINATION lib)

# Tests
add_subdirectory(tests)
