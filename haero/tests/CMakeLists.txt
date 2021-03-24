include(EkatCreateUnitTest)
link_directories(${EKAT_LIBRARY_DIR}) # Needed to link unit tests.

# This simple unit test gets us off the ground.
EkatCreateUnitTest(hello_world hello_world.cpp
                   LIBS ${HAERO_LIBRARIES})

#----------------------------------------------------
# Unit tests for basic data structures and functions
#----------------------------------------------------

EkatCreateUnitTest(view_pack_helpers view_helper_tests.cpp
                   LIBS ${HAERO_LIBRARIES})

EkatCreateUnitTest(floating_point_tests floating_point_tests.cpp
                   LIBS ${HAERO_LIBRARIES})

EkatCreateUnitTest(mode_tests mode_tests.cpp
                   LIBS ${HAERO_LIBRARIES})

EkatCreateUnitTest(species_tests species_tests.cpp
                   LIBS ${HAERO_LIBRARIES})

EkatCreateUnitTest(model_tests model_tests.cpp
                   LIBS ${HAERO_LIBRARIES})

EkatCreateUnitTest(prognostics_tests prognostics_tests.cpp
                   LIBS ${HAERO_LIBRARIES})

EkatCreateUnitTest(tendencies_tests tendencies_tests.cpp
                   LIBS ${HAERO_LIBRARIES})

add_library(chemUtil chemUtil.cpp)
add_dependencies(chemUtil tchem)
EkatCreateUnitTest(tchem_batch_reactor tchem_batch_reactor.cpp
                   LIBS chemUtil;${HAERO_LIBRARIES})

EkatCreateUnitTest(kokkos_team_execute kokkos_team_execute.cpp
                   LIBS ${HAERO_LIBRARIES})

EkatCreateUnitTest(process_tests process_tests.cpp
                   LIBS ${HAERO_LIBRARIES})

#----------------------------------------------------
# Unit tests for prognostic processes
#----------------------------------------------------

if (HAERO_FORTRAN)
  # Fortran process stub tests -- checks that data can be passed between C++ and
  # Fortran
  add_library(haero_fprocess_stubs
              prog_fprocess_stub.cpp prog_fprocess_stub.F90
              diag_fprocess_stub.cpp diag_fprocess_stub.F90
              mam_nucleation_test_bridge.F90)
  add_dependencies(haero_fprocess_stubs haero)
  EkatCreateUnitTest(fprocess_stub_tests fprocess_stub_tests.cpp
                     LIBS haero_fprocess_stubs;${HAERO_LIBRARIES})

  EkatCreateUnitTest(mam_nucleation_process_tests mam_nucleation_process_tests.cpp
                     LIBS haero_fprocess_stubs;${HAERO_LIBRARIES})
endif()

#---------------------------
# Aerosol process unit tests
#---------------------------

if (HAERO_FORTRAN)
  EkatCreateUnitTest(mam_nucleation_fprocess_tests mam_nucleation_fprocess_tests.cpp
                     LIBS ${HAERO_LIBRARIES})
endif()
