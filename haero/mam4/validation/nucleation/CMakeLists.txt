set(NUC_VALIDATION_DIR ${MAM_X_VALIDATION_DIR}/nucleation)
set(NUC_BASELINE_DIR ${MAM4BOX_BASELINE_DIR}/nucleation)

# These subdirectories contain Skywalker drivers for MAM4 parameterizations.
# Include directory for .mod files.
include_directories(${PROJECT_BINARY_DIR}/haero)
include_directories(${PROJECT_BINARY_DIR}/validation)

# We use a single driver for all nucleation-related parameterizations.
add_executable(nucleation_driver nucleation_driver.cpp
               nucleation_rate.cpp
               nucleation_thresh.cpp)
target_link_libraries(nucleation_driver skywalker;validation;${HAERO_LIBRARIES})

# Run the drivers in several configurations to produce datasets.

# Nucleation rate
foreach (input vehkamaki2002_contour;vehkamaki2002_fig8;
               vehkamaki2002_fig9;vehkamaki2002_fig10;vehkamaki2002_fig11;
               merikanto2007_fig2;merikanto2007_fig3;merikanto2007_fig4)
  add_test(validate_mam4_nucleation_rate_${input} nucleation_driver ${NUC_VALIDATION_DIR}/${input}.yaml)
endforeach()

# Nucleation threshold
foreach (input vehkamaki2002_fig7)
  add_test(validate_mam4_nucleation_thresh_${input} nucleation_driver ${NUC_VALIDATION_DIR}/${input}.yaml)
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
