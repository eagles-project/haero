add_library(toy_problem toy_problem.cpp)

add_dependencies(toy_problem chem_standalone)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/toy_problem_test.cpp.in toy_problem_test.cpp)
EkatCreateUnitTest(toy_problem_test ${CMAKE_CURRENT_BINARY_DIR}/toy_problem_test.cpp LIBS chem_standalone toy_problem ${HAERO_LIBRARIES} ${HAERO_DRIVER_LIBRARIES})
