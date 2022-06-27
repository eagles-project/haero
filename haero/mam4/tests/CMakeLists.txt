include(create_unit_test)
#link_directories(${EKAT_LIBRARY_DIR}) # Needed to link unit tests.

CreateUnitTest(mam4_nucleation mam4_nucleation_tests.cpp)

