set(NUC_VALIDATION_DIR ${MAM_X_VALIDATION_DIR}/nucleation)

# Include directory for .mod files.
include_directories(${PROJECT_BINARY_DIR}/haero)
include_directories(${PROJECT_BINARY_DIR}/validation)

# These skywalker drivers compute binary and ternary nucleation rates.
add_executable(nucleation_rate nucleation_rate.cpp)
target_link_libraries(nucleation_rate skywalker;validation;${HAERO_LIBRARIES})

# This driver computes the threshold concentration of H2SO4 at given T and RH
# for nucleation to occur.
add_executable(nucleation_thresh nucleation_thresh.cpp)
target_link_libraries(nucleation_thresh skywalker;validation;${HAERO_LIBRARIES})

# This driver computes the growth SO4 and NH4 aerosol particeles from
# H2SO4 and NH3 gas at given T and RH.
add_executable(nucleation_growth nucleation_growth.cpp)
target_link_libraries(nucleation_growth skywalker;validation;${HAERO_LIBRARIES})

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

# Nucleation growth
#foreach (input nuc_growth_case1 nuc_growth_case2 nuc_growth_case3)
#  add_test(gen_${input} nucleation_thresh ${NUC_VALIDATION_DIR}/${input}.yaml)
#endforeach()

# Copy generate_plots.py to our binary directory for ease of use.
configure_file(
  ${NUC_VALIDATION_DIR}/generate_plots.py
  ${CMAKE_CURRENT_BINARY_DIR}/generate_plots.py
  COPYONLY
)
