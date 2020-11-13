include(EkatCreateUnitTest)

# This simple unit test gets us off the ground.
EkatCreateUnitTest(hello_world hello_world.cpp
                   LIBS ${HAERO_LIBRARIES})

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

EkatCreateUnitTest(prog_process_stub_tests prog_process_stub_tests.cpp
                   LIBS ${HAERO_LIBRARIES})

EkatCreateUnitTest(diag_process_stub_tests diag_process_stub_tests.cpp
                   LIBS ${HAERO_LIBRARIES})
