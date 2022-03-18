include(EkatCreateUnitTest)
link_directories(${EKAT_LIBRARY_DIR}) # Needed to link unit tests.

EkatCreateUnitTest(chem_driver_tests chem_driver_tests.cpp LIBS
                                                           chem_standalone
                                                           yaml_cpp
                                                           ${HAERO_LIBRARIES})
