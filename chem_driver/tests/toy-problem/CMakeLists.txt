add_library(toy_problem toy_problem.cpp)

add_dependencies(toy_problem chem_standalone)

include(EkatCreateUnitTest)
EkatCreateUnitTest(toy_problem_test toy_problem_test.cpp
                   LIBS toy_problem;chem_standalone;yaml_cpp;${HAERO_LIBRARIES})
