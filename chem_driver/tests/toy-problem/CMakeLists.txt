include(create_unit_test)

add_library(toy_problem toy_problem.cpp)

add_dependencies(toy_problem chem_standalone)

CreateUnitTest(toy_problem_test toy_problem_test.cpp LIBS chem_standalone toy_problem ${HAERO_LIBRARIES} ${HAERO_DRIVER_LIBRARIES})
