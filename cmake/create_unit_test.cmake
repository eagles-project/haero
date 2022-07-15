# Use this function to create simple HAERO unit tests (serial, nothing fancy).
# Use EkatCreateUnitTest if you need something more complicated.

include(EkatCreateUnitTest)

function(CreateUnitTest test_name test_srcs)
  if (HAERO_MPI_EXEC)
    if (HAERO_MPI_EXTRA_ARGS)
      EkatCreateUnitTest(${test_name} ${test_srcs}
        MPI_EXEC_NAME ${HAERO_MPI_EXEC}
        MPI_NP_FLAG ${HAERO_MPI_NP_FLAG}
        MPI_EXTRA_ARGS ${HAERO_MPI_EXTRA_ARGS}
        LIBS ${HAERO_LIBRARIES})
    else()
      EkatCreateUnitTest(${test_name} ${test_srcs}
        MPI_EXEC_NAME ${HAERO_MPI_EXEC}
        MPI_NP_FLAG ${HAERO_MPI_NP_FLAG}
        LIBS ${HAERO_LIBRARIES})
    endif()
  else()
    EkatCreateUnitTest(${test_name} ${test_srcs} LIBS ${HAERO_LIBRARIES})
  endif()
endfunction()

function(CreateSkywalkerTest test_name test_srcs)
  if (HAERO_MPI_EXEC)
    if (HAERO_MPI_EXTRA_ARGS)
      EkatCreateUnitTest(${test_name} ${test_srcs}
        MPI_EXEC_NAME ${HAERO_MPI_EXEC}
        MPI_NP_FLAG ${HAERO_MPI_NP_FLAG}
        MPI_EXTRA_ARGS ${HAERO_MPI_EXTRA_ARGS}
        LIBS skywalker;${HAERO_LIBRARIES})
    else()
      EkatCreateUnitTest(${test_name} ${test_srcs}
        MPI_EXEC_NAME ${HAERO_MPI_EXEC}
        MPI_NP_FLAG ${HAERO_MPI_NP_FLAG}
        LIBS skywalker;${HAERO_LIBRARIES})
    endif()
  else()
    EkatCreateUnitTest(${test_name} ${test_srcs} LIBS skywalker;${HAERO_LIBRARIES})
  endif()
endfunction()
