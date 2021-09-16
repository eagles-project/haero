# Build the validation program.
add_executable(validate validate.cpp)
target_link_libraries(validate skywalker;yaml_cpp;${HAERO_LIBRARIES})

# Copy generate_plots.py to our binary directory for ease of use.
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/generate_plots.py
  ${CMAKE_CURRENT_BINARY_DIR}/generate_plots.py
  COPYONLY
)

# Run it in several configurations to produce datasets.
add_test(validate_binary_nuc_rate validate ${CMAKE_CURRENT_SOURCE_DIR}/binary_nuc_contour.yaml)
set_tests_properties(validate_binary_nuc_rate PROPERTIES
                     WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
