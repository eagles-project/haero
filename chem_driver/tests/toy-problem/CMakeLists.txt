add_library(toy_problem toy_problem.cpp)

add_dependencies(toy_problem chem_standalone)

EkatCreateUnitTest(toy_problem_test toy_problem_test.cpp LIBS chem_standalone toy_problem ${HAERO_LIBRARIES} ${HAERO_DRIVER_LIBRARIES})

# install(FILES toy_problem_input.yml DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")
