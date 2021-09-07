if (HAERO_FORTRAN)

  # Make sure Skywalker's Fortran 90 interface behaves properly.
  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/skywalker_f90_test.F90.in
    ${CMAKE_CURRENT_BINARY_DIR}/skywalker_f90_test.F90
    @ONLY
  )
  add_executable(skywalker_f90_test ${CMAKE_CURRENT_BINARY_DIR}/skywalker_f90_test.F90)
  target_link_libraries(skywalker_f90_test skywalker;yaml_cpp;${HAERO_LIBRARIES})
  add_test(skywalker_f90_test ${CMAKE_CURRENT_BINARY_DIR}/skywalker_f90_test)
  set_tests_properties(skywalker_f90_test PROPERTIES
    PASS_REGULAR_EXPRESSION "PASS")

  # This tests a user-defined parameter walk (no atmosphere/aerosol/gas input).
  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/skywalker_f90_user_test.F90.in
    ${CMAKE_CURRENT_BINARY_DIR}/skywalker_f90_user_test.F90
    @ONLY
  )
  add_executable(skywalker_f90_user_test ${CMAKE_CURRENT_BINARY_DIR}/skywalker_f90_user_test.F90)
  target_link_libraries(skywalker_f90_user_test skywalker;yaml_cpp;${HAERO_LIBRARIES})
  add_test(skywalker_f90_user_test ${CMAKE_CURRENT_BINARY_DIR}/skywalker_f90_user_test)
  set_tests_properties(skywalker_f90_user_test PROPERTIES
    PASS_REGULAR_EXPRESSION "PASS")
endif()
