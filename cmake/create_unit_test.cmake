# Use this function to create simple HAERO unit tests (serial, nothing fancy).
# Use EkatCreateUnitTest if you need something more complicated.

include(EkatCreateUnitTest)

function(CreateUnitTest test_name test_srcs)
  # Some machines always need mpiexec (or srun, perhaps), even for nproc = 1.
  # TODO: Use srun where needed (for this we must detect SLURM).
  EkatCreateUnitTest(${test_name} ${test_srcs}
                     MPI_EXEC_NAME mpiexec
                     MPI_NP_FLAG -n
                     LIBS ${HAERO_LIBRARIES})
endfunction()
