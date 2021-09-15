include(create_unit_test)
link_directories(${EKAT_LIBRARY_DIR}) # Needed to link unit tests.

set(CREATE_ATM_CONDITIONS_TEST "
  CreateUnitTest(atm_conditions atm_condition_test.cpp LIBS haero_standalone ${HAERO_LIBRARIES}
   ${HAERO_DRIVER_LIBRARIES} PROPERTIES PASS_REGULAR_EXPRESSION \"zp1000 = 1000\"
  ")
if (HAERO_MPI_EXEC)
  string(APPEND CREATE_ATM_CONDITIONS_TEST "
    MPI_EXEC_NAME ${HAERO_MPI_EXEC}
    MPI_NP_FLAG ${HAERO_MPI_NP_FLAG}
    ")
endif()

if (HAERO_MPI_EXTRA_ARGS)
  string(APPEND CREATE_ATM_CONDITIONS_TEST "
    MPI_EXTRA_ARGS ${HAERO_MPI_EXTRA_ARGS}
    ")
else()

string(APPEND CREATE_ATM_CONDITIONS_TEST ")")

cmake_language(EVAL CODE "${CREATE_ATM_CONDITIONS_TEST}")

# YAML input tests.
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/read_yaml_test.cpp.in read_yaml_test.cpp)
CreateUnitTest(read_yaml ${CMAKE_CURRENT_BINARY_DIR}/read_yaml_test.cpp LIBS haero_standalone ${HAERO_LIBRARIES}
  ${HAERO_DRIVER_LIBRARIES})

CreateUnitTest(ncwriter ncwriter_tests.cpp LIBS haero_standalone ${HAERO_LIBRARIES} ${HAERO_DRIVER_LIBRARIES})

CreateUnitTest(host_dynamics host_dynamics_test.cpp LIBS haero_standalone ${HAERO_LIBRARIES}
  ${HAERO_DRIVER_LIBRARIES})

#CreateUnitTest(host_model host_model_tests.cpp LIBS haero_standalone ${HAERO_LIBRARIES} ${HAERO_DRIVER_LIBRARIES})

#CreateUnitTest(host_dynamics_tendencies host_dynamics_tendency_tests.cpp LIBS haero_standalone
#  ${HAERO_DRIVER_LIBRARIES} ${HAERO_LIBRARIES})
