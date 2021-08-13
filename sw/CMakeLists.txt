include_directories(${CMAKE_CURRENT_BINARY_DIR} ${HAERO_INCDIRS} ${HAERO_DRIVER_INCDIRS})
include_directories(${PROJECT_BINARY_DIR}/haero)

# Generate skywalker.hpp with the correct floating point precision type.
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/skywalker.hpp.in
  ${CMAKE_CURRENT_BINARY_DIR}/skywalker.hpp
  @ONLY
)

# Here we build a skywalker library that can be used by our MAM box model.
if (HAERO_FORTRAN)
  set(SW_FORTRAN_SOURCES skywalker.F90)
endif()
add_library(skywalker skywalker.cpp load_ensemble.cpp write_py_module.cpp
            ${SW_FORTRAN_SOURCES})
set_target_properties(skywalker PROPERTIES OUTPUT_NAME skywalker_${HAERO_REAL_TYPE})
add_dependencies(skywalker haero)

# Here's our Haero skywalker executable.
add_executable(skywalker_exe main.cpp)
set_target_properties(skywalker_exe PROPERTIES OUTPUT_NAME skywalker)
target_link_libraries(skywalker_exe skywalker;yaml_cpp;${HAERO_LIBRARIES})

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/skywalker.hpp
        DESTINATION include/skywalker)
if (HAERO_FORTRAN)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/skywalker.mod
                ${PROJECT_BINARY_DIR}/haero/haero.mod
                ${PROJECT_BINARY_DIR}/haero/haero_precision.mod
                ${PROJECT_BINARY_DIR}/haero/haero_constants.mod
          DESTINATION include/skywalker)
endif()
install(TARGETS skywalker DESTINATION lib)
install(TARGETS skywalker_exe DESTINATION bin)

# Generate skywalker.cmake with installation information
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/skywalker.cmake.in
  ${PROJECT_BINARY_DIR}/cmake/skywalker.cmake
  @ONLY
)
install(FILES ${PROJECT_BINARY_DIR}/cmake/skywalker.cmake
        DESTINATION cmake)

# Skywalker smoke test -- to ensure that we can build and run skywalker!
add_test(skywalker_smoke_test ${CMAKE_CURRENT_BINARY_DIR}/skywalker ${CMAKE_CURRENT_SOURCE_DIR}/nucleation.yaml)
set_tests_properties(skywalker_smoke_test PROPERTIES
  FAIL_REGULAR_EXPRESSION "error:"
)
