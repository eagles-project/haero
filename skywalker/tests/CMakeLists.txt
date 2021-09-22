if (HAERO_FORTRAN)

  # Make sure Skywalker's Fortran 90 interface behaves properly.
  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/f90_test.yaml
    ${CMAKE_CURRENT_BINARY_DIR}/f90_test.yaml
    COPYONLY
  )
  add_executable(skywalker_f90_test skywalker_f90_test.F90)
  target_link_libraries(skywalker_f90_test skywalker;yaml_cpp;${HAERO_LIBRARIES})
  set_target_properties(skywalker_f90_test PROPERTIES
    LINKER_LANGUAGE Fortran)
  add_test(skywalker_f90_test skywalker_f90_test)
  set_tests_properties(skywalker_f90_test PROPERTIES
    PASS_REGULAR_EXPRESSION "PASS")

  # This tests a user-defined parameter walk (no atmosphere/aerosol/gas input).
  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/f90_user_test.yaml
    ${CMAKE_CURRENT_BINARY_DIR}/f90_user_test.yaml
    COPYONLY
  )
  add_executable(skywalker_f90_user_test skywalker_f90_user_test.F90)
  target_link_libraries(skywalker_f90_user_test skywalker;yaml_cpp;${HAERO_LIBRARIES})
  set_target_properties(skywalker_f90_user_test PROPERTIES
    LINKER_LANGUAGE Fortran)
  add_test(skywalker_f90_user_test skywalker_f90_user_test)
  set_tests_properties(skywalker_f90_user_test PROPERTIES
    PASS_REGULAR_EXPRESSION "PASS")
endif()
