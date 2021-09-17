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
foreach (input vehkamaki2002_contour;vehkamaki2002_fig7;vehkamaki2002_fig8;
               vehkamaki2002_fig9;vehkamaki2002_fig10;vehkamaki2002_fig11;
               merikanto2007_fig2;merikanto2007_fig3;merikanto2007_fig4)
  add_test(gen_${input} validate ${CMAKE_CURRENT_SOURCE_DIR}/${input}.yaml)
endforeach()
