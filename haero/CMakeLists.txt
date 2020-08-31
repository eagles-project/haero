include_directories(${CMAKE_CURRENT_BINARY_DIR})

# Generate a file containing our config parameters.
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/haero_config.hpp.in
  ${CMAKE_CURRENT_BINARY_DIR}/haero_config.hpp
  @ONLY
)

# Generate haero_version.cpp with proper version/revision information.
add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/haero_version.cpp
  COMMAND ${BASH} ${PROJECT_SOURCE_DIR}/tools/update_version_info.sh ${PROJECT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/haero_version.cpp
  COMMENT "Generating haero version information..."
)

# Update the haero version information if it's already been generated
add_custom_target(update_version_info ALL
  COMMAND ${BASH} ${PROJECT_SOURCE_DIR}/tools/update_version_info.sh ${PROJECT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/haero_version.cpp
  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/haero_version.cpp
  COMMENT "Updating haero version information..."
)

add_library(haero
            ${CMAKE_CURRENT_BINARY_DIR}/haero_version.cpp
            )
set(HAERO_LIBRARIES haero;${HAERO_LIBRARIES})
set(HAERO_LIBRARIES haero;${HAERO_LIBRARIES} PARENT_SCOPE)

# Installation targets
install(TARGETS haero DESTINATION lib)

# Tests
#add_subdirectory(tests)

