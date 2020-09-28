include_directories(${CMAKE_CURRENT_BINARY_DIR})

# Generate a file containing our config parameters.
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

add_library(haero
            ${CMAKE_CURRENT_BINARY_DIR}/haero_version.cpp
            utils.cpp
            view_pack_helpers.cpp
            )
add_dependencies(haero ekat update_version_info)
set(HAERO_LIBRARIES haero;${HAERO_LIBRARIES})
set(HAERO_LIBRARIES haero;${HAERO_LIBRARIES} PARENT_SCOPE)

# Installation targets
install(TARGETS haero DESTINATION lib)
install(FILES haero.hpp
              ${CMAKE_CURRENT_BINARY_DIR}/haero_config.hpp
              mode.hpp
              physical_constants.hpp
              utils.hpp
              view_pack_helpers.hpp
              species.hpp
              utils.hpp
              floating_point.hpp
        DESTINATION include/haero)

# Library unit tests.
add_subdirectory(tests)

