set(NUC_VALIDATION_DIR ${MAM_X_VALIDATION_DIR}/gasaerexch)

# Include directory for .mod files.
include_directories(${PROJECT_BINARY_DIR}/haero)
include_directories(${PROJECT_BINARY_DIR}/validation)


# These skywalker drivers compute gasaerexch update rates.
add_executable(exe_skywalker_gasaerexch_uptkrates_1box1gas  skywalker_gasaerexch_uptkrates_1box1gas.F90)
target_link_libraries(exe_skywalker_gasaerexch_uptkrates_1box1gas skywalker_f90;skywalker;validation_f90;${HAERO_LIBRARIES})

# Run the drivers in several configurations to produce datasets.

# gasaerexch rate
foreach (input skywalker_gasaerexch_uptkrates_1box1gas)
  add_test(gen_${input} exe_skywalker_gasaerexch_uptkrates_1box1gas ${NUC_VALIDATION_DIR}/${input}.yaml)
endforeach()

# gasaerexch growth
#foreach (input nuc_growth_case1 nuc_growth_case2 nuc_growth_case3)
#  add_test(gen_${input} gasaerexch_thresh ${NUC_VALIDATION_DIR}/${input}.yaml)
#endforeach()

# Copy generate_plots.py to our binary directory for ease of use.
# configure_file(
#   ${NUC_VALIDATION_DIR}/generate_plots.py
#   ${CMAKE_CURRENT_BINARY_DIR}/generate_plots.py
#   COPYONLY
# )
