# Generates haero.cmake (for downstream aerosol packages) and install it.
macro(HaeroGenerateConfig)
  # generate the cmake file used by mam4xx
  configure_file(
    ${PROJECT_SOURCE_DIR}/cmake/haero.cmake.in
    ${PROJECT_BINARY_DIR}/share/haero.cmake
    @ONLY
  )
  install(DIRECTORY ${PROJECT_BINARY_DIR}/share DESTINATION .)
endmacro()
