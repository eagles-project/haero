include(EkatCreateUnitTest)

EkatCreateUnitTest(phys_column phys_column_test.cpp LIBS haero_standalone ${HAERO_LIBRARIES})
EkatCreateUnitTest(atm_conditons atm_condition_test.cpp LIBS haero_standalone ${HAERO_LIBRARIES})
EkatCreateUnitTest(dynamics_column dyn_column_test.cpp LIBS haero_standalone ${HAERO_LIBRARIES})
