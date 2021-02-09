include_directories(${CMAKE_CURRENT_BINARY_DIR})

# Generate files containing our config parameters.
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/haero_config.hpp.in
  ${CMAKE_CURRENT_BINARY_DIR}/haero_config.hpp
  @ONLY
)

# Generate a blank haero_version.cpp file to get us off the ground.
add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/haero_version.cpp
  COMMAND ${CMAKE_COMMAND} -E touch ${CMAKE_CURRENT_BINARY_DIR}/haero_version.cpp
  COMMENT "Generating initial haero_version.cpp..."
)

# Update the haero version information if it's already been generated
add_custom_target(update_version_info ALL
  COMMAND ${BASH} ${PROJECT_SOURCE_DIR}/tools/update_version_info.sh ${PROJECT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/haero_version.cpp
  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/haero_version.cpp
  COMMENT "Updating haero version information..."
)

# Gather Haero process implementations
add_subdirectory(processes)

if (HAERO_FORTRAN)
  set(HAERO_FORTRAN_SOURCES haero.F90)
endif()

add_library(haero
            ${CMAKE_CURRENT_BINARY_DIR}/haero_version.cpp
            ${HAERO_FORTRAN_SOURCES}
            atmosphere.cpp
            chemistry.cpp
            diagnostics.cpp
            model.cpp
            prognostics.cpp
            selected_processes.cpp
            tendencies.cpp
            utils.cpp
            view_pack_helpers.cpp
            ${HAERO_PROCESSES} # populated by the processes/ ѕubdirectory
            )
add_dependencies(haero ekat update_version_info tchem)
set(HAERO_LIBRARIES haero;${HAERO_LIBRARIES})
set(HAERO_LIBRARIES haero;${HAERO_LIBRARIES} PARENT_SCOPE)

# Installation targets
install(TARGETS haero DESTINATION lib)
install(FILES haero.hpp
              ${CMAKE_CURRENT_BINARY_DIR}/haero_config.hpp
              diagnostics.hpp
              floating_point.hpp
              mode.hpp
              model.hpp
              selected_processes.hpp
              physical_constants.hpp
              prognostics.hpp
              species.hpp
              utils.hpp
              view_pack_helpers.hpp
        DESTINATION include/haero)

# Library unit tests.
add_subdirectory(tests)

