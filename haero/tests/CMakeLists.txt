include(EkatCreateUnitTest)

#------------
# Unit Tests
#------------

# We use the EkatCreateUnitTest CMake function to create unit tests that are
# configured properly to work in the Kokkos environment provided by EKAT.
# EkatCreateUnitTest must be called with the following options:
# 1. LIBS ${HAERO_LIBRARIES} <-- links against Haero
# 2. EXCLUDE_TEST_SESSION    <-- uses Haero's setup/breakdown functions

EkatCreateUnitTest(math_tests math_tests.cpp
                   LIBS ${HAERO_LIBRARIES} EXCLUDE_TEST_SESSION)
EkatCreateUnitTest(testing_tests testing_tests.cpp
                   LIBS ${HAERO_LIBRARIES} EXCLUDE_TEST_SESSION)
