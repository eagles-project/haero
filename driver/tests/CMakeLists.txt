include(EkatCreateUnitTest)

EkatCreateUnitTest(atm_conditions atm_condition_test.cpp LIBS haero_standalone ${HAERO_LIBRARIES}
 ${HAERO_DRIVER_LIBRARIES} PROPERTIES PASS_REGULAR_EXPRESSION "zp1000 = 1000")

# YAML input tests.
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/read_yaml_test.cpp.in read_yaml_test.cpp)
EkatCreateUnitTest(read_yaml ${CMAKE_CURRENT_BINARY_DIR}/read_yaml_test.cpp LIBS haero_standalone ${HAERO_LIBRARIES}  
  ${HAERO_DRIVER_LIBRARIES})

EkatCreateUnitTest(ncwriter ncwriter_tests.cpp LIBS haero_standalone ${HAERO_LIBRARIES} ${HAERO_DRIVER_LIBRARIES})

EkatCreateUnitTest(host_dynamics host_dynamics_test.cpp LIBS haero_standalone ${HAERO_LIBRARIES}
  ${HAERO_DRIVER_LIBRARIES})

#EkatCreateUnitTest(dynamics_column dyn_column_test.cpp LIBS haero_standalone ${HAERO_DRIVER_LIBRARIES} ${HAERO_LIBRARIES} )
#EkatCreateUnitTest(phys_column phys_column_test.cpp LIBS haero_standalone ${HAERO_LIBRARIES})
