include(create_unit_test)
#link_directories(${EKAT_LIBRARY_DIR}) # Needed to link unit tests.

CreateUnitTest(mam4_nucleation_unit_tests mam4_nucleation_unit_tests.cpp)

