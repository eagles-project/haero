include(EkatCreateUnitTest)

EkatCreateUnitTest(math_tests math_tests.cpp
                   LIBS ${HAERO_LIBRARIES} EXCLUDE_TEST_SESSION)
EkatCreateUnitTest(testing_tests testing_tests.cpp
                   LIBS ${HAERO_LIBRARIES} EXCLUDE_TEST_SESSION)
