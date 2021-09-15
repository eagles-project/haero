# Use this function to create simple HAERO unit tests (serial, nothing fancy).
# Use EkatCreateUnitTest if you need something more complicated.

include(EkatCreateUnitTest)

function(CreateUnitTest test_name test_srcs)

  # Only parse arguments after the positional args
  cmake_parse_arguments(CreateUnitTest "" "" "LIBS" ${ARGN})

  set(CreateUnitTest_LIBRARIES ${HAERO_LIBRARIES} ${CreateUnitTest_LIBS})

  if (HAERO_MPI_EXEC)
    if (HAERO_MPI_EXTRA_ARGS)
      EkatCreateUnitTest(${test_name} ${test_srcs}
        MPI_EXEC_NAME ${HAERO_MPI_EXEC}
        MPI_NP_FLAG ${HAERO_MPI_NP_FLAG}
        MPI_EXTRA_ARGS ${HAERO_MPI_EXTRA_ARGS}
        LIBS ${CreateUnitTest_LIBRARIES})
    else()
      EkatCreateUnitTest(${test_name} ${test_srcs}
        MPI_EXEC_NAME ${HAERO_MPI_EXEC}
        MPI_NP_FLAG ${HAERO_MPI_NP_FLAG}
        LIBS ${CreateUnitTest_LIBRARIES})
    endif()
  else()
    EkatCreateUnitTest(${test_name} ${test_srcs} LIBS ${CreateUnitTest_LIBRARIES})
  endif()
endfunction()
