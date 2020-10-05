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

EkatCreateUnitTest(aero_model_tests aero_model_tests.cpp
                   LIBS ${HAERO_LIBRARIES})

EkatCreateUnitTest(aero_state_tests aero_state_tests.cpp
                   LIBS ${HAERO_LIBRARIES})

EkatCreateUnitTest(aero_tendencies_tests aero_tendencies_tests.cpp
                   LIBS ${HAERO_LIBRARIES})
