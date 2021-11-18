include(EkatCreateUnitTest)
link_directories(${EKAT_LIBRARY_DIR}) # Needed to link unit tests.

add_subdirectory(util)

EkatCreateUnitTest(chem_driver_tests chem_driver_tests.cpp LIBS chem_standalone read_solver_params ${HAERO_LIBRARIES})
