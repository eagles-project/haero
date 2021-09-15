include(create_unit_test)
link_directories(${EKAT_LIBRARY_DIR}) # Needed to link unit tests.

CreateUnitTest(atm_conditions atm_condition_test.cpp LIBS haero_standalone ${HAERO_LIBRARIES}
 ${HAERO_DRIVER_LIBRARIES} PROPERTIES PASS_REGULAR_EXPRESSION "zp1000 = 1000")

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
