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

EkatCreateUnitTest(aerosol_species_tests aerosol_species_tests.cpp
                   LIBS ${HAERO_LIBRARIES})

EkatCreateUnitTest(gas_species_tests gas_species_tests.cpp
                   LIBS ${HAERO_LIBRARIES})

EkatCreateUnitTest(model_tests model_tests.cpp
                   LIBS ${HAERO_LIBRARIES})

EkatCreateUnitTest(prognostics_tests prognostics_tests.cpp
                   LIBS ${HAERO_LIBRARIES})

EkatCreateUnitTest(tendencies_tests tendencies_tests.cpp
                   LIBS ${HAERO_LIBRARIES})

EkatCreateUnitTest(kokkos_team_execute kokkos_team_execute.cpp
                   LIBS ${HAERO_LIBRARIES})

EkatCreateUnitTest(process_tests process_tests.cpp
                   LIBS ${HAERO_LIBRARIES})

EkatCreateUnitTest(cuda_tests cuda_tests.cpp
                   LIBS ${HAERO_LIBRARIES})

EkatCreateUnitTest(math_tests math_helper_tests.cpp
                   LIBS ${HAERO_LIBRARIES})

EkatCreateUnitTest(kohler_tests kohler_tests.cpp
                   LIBS ${HAERO_LIBRARIES})

#EkatCreateUnitTest(modal_wet_radius modal_wet_radius_tests.cpp
#                   LIBS ${HAERO_LIBRARIES})

#----------------------------------------------------
# Unit tests for aerosol processes
#----------------------------------------------------

EkatCreateUnitTest(mam_nucleation_process_tests mam_nucleation_process_tests.cpp
                   LIBS ${HAERO_LIBRARIES})

EkatCreateUnitTest(mam_nucleation_packed_tests mam_nucleation_packed_tests.cpp
                   LIBS ${HAERO_LIBRARIES})

if (HAERO_FORTRAN)
  # Fortran process stub tests -- checks that data can be passed between C++ and
  # Fortran
  add_library(haero_faerosol_process_stubs
              faerosol_process_stub.cpp faerosol_process_stub.F90
              mam_nucleation_test_bridge.F90)
  add_dependencies(haero_faerosol_process_stubs haero)
  EkatCreateUnitTest(faerosol_process_stub_tests faerosol_process_stub_tests.cpp
                     LIBS haero_faerosol_process_stubs;${HAERO_LIBRARIES})

  # TODO: the Fortran version of the nucleation process doesn't currently pass
  # TODO: tests in single-precision builds. We'll punt for now, but should
  # TODO: probably try to figure this out soon.
  if (HAERO_DOUBLE_PRECISION)
    # Aerosol process unit tests
    EkatCreateUnitTest(mam_nucleation_fprocess_tests mam_nucleation_fprocess_tests.cpp
                       LIBS haero_faerosol_process_stubs;${HAERO_LIBRARIES})
  endif()
endif()

#----------------------------------------------------
# Chemistry unit tests
#----------------------------------------------------

if (HAERO_DOUBLE_PRECISION)
  add_library(chemUtil chemUtil.cpp)
  add_dependencies(chemUtil tchem)
  EkatCreateUnitTest(tchem_batch_reactor_tests tchem_batch_reactor_tests.cpp
                     LIBS chemUtil;${HAERO_LIBRARIES})
endif()
