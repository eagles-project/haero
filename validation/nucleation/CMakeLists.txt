set(NUC_VALIDATION_DIR ${MAM_X_VALIDATION_DIR}/nucleation)

# This skywalker driver computes binary and ternary nucleation rates.
add_executable(nucleation_rate nucleation_rate.cpp)
target_link_libraries(nucleation_rate skywalker;${HAERO_LIBRARIES})

# This driver computes the threshold concentration of H2SO4 at given T and RH
# for nucleation to occur.
add_executable(nucleation_thresh nucleation_thresh.cpp)
target_link_libraries(nucleation_thresh skywalker;${HAERO_LIBRARIES})

# Run the drivers in several configurations to produce datasets.

# Nucleation rate
foreach (input vehkamaki2002_contour;vehkamaki2002_fig8;
               vehkamaki2002_fig9;vehkamaki2002_fig10;vehkamaki2002_fig11;
               merikanto2007_fig2;merikanto2007_fig3;merikanto2007_fig4)
  add_test(gen_${input} nucleation_rate ${NUC_VALIDATION_DIR}/${input}.yaml)
endforeach()

# Nucleation threshold
foreach (input vehkamaki2002_fig7)
  add_test(gen_${input} nucleation_thresh ${NUC_VALIDATION_DIR}/${input}.yaml)
endforeach()

# Copy generate_plots.py to our binary directory for ease of use.
configure_file(
  ${NUC_VALIDATION_DIR}/generate_plots.py
  ${CMAKE_CURRENT_BINARY_DIR}/generate_plots.py
  COPYONLY
)
