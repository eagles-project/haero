include(EkatCreateUnitTest)
link_directories(${EKAT_LIBRARY_DIR}) # Needed to link unit tests.

EkatCreateUnitTest(atm_conditions atm_condition_test.cpp LIBS haero_standalone ${HAERO_LIBRARIES})

# YAML input tests.
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/read_yaml_test.cpp.in read_yaml_test.cpp)
EkatCreateUnitTest(read_yaml read_yaml_test.cpp LIBS haero_standalone ${HAERO_LIBRARIES})

EkatCreateUnitTest(ncwriter ncwriter_tests.cpp LIBS haero_standalone ${HAERO_LIBRARIES} ${HAERO_DRIVER_LIBRARIES})
EkatCreateUnitTest(dynamics_column dyn_column_test.cpp LIBS haero_standalone ${HAERO_DRIVER_LIBRARIES} ${HAERO_LIBRARIES} )
#EkatCreateUnitTest(phys_column phys_column_test.cpp LIBS haero_standalone ${HAERO_LIBRARIES})
