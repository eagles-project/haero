set(VALIDATION_DIR ${MAM_X_VALIDATION_DIR}/gasaerexch)

# Include directory for .mod files.
include_directories(${PROJECT_BINARY_DIR}/haero)
include_directories(${PROJECT_BINARY_DIR}/validation)


# These skywalker drivers compute gasaerexch update rates.
add_executable(exe_skywalker_gasaerexch_uptkrates_1box1gas  skywalker_gasaerexch_uptkrates_1box1gas.F90)
target_link_libraries(exe_skywalker_gasaerexch_uptkrates_1box1gas skywalker_f90;skywalker;validation_f90;${HAERO_LIBRARIES})

# Run the drivers in several configurations to produce datasets.

# gasaerexch rate
foreach (input skywalker_gasaerexch_1subarea_1gas_nonvolatile
               skywalker_gasaerexch_uptkrates_1box1gas
               skywalker_gasaerexch_uptkrates_1box1gas_plot
               skywalker_gasaerexch_veh02_fig08)
  add_test(gen_${input} exe_skywalker_gasaerexch_uptkrates_1box1gas ${VALIDATION_DIR}/${input}.yaml)
endforeach()
