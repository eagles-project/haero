
# Design document.
find_package(LATEX COMPONENTS PDFLATEX BIBTEX)
if (LATEX_PDFLATEX_FOUND AND LATEX_BIBTEX_FOUND)
  add_subdirectory(design)
  #
  # Putting the .tex configure file in both source and binary dirs
  # to allow non-CMake compilation of the .tex docs in-source
  #
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/haero_design_inc.tex.in
                 ${CMAKE_CURRENT_BINARY_DIR}/haero_design_inc.tex
                 @ONLY
  )
 configure_file(${CMAKE_CURRENT_SOURCE_DIR}/haero_design_inc.tex.in
                ${CMAKE_CURRENT_SOURCE_DIR}/haero_design_inc.tex
                @ONLY
  )
endif()

# Doxygen reference documentation
find_package(Doxygen OPTIONAL_COMPONENTS dot QUIET)
if (DOXYGEN_FOUND)
  add_subdirectory(doxygen)
endif()

# Create a "docs" target that builds all documentation.
add_custom_target(docs DEPENDS ${DOC_TARGETS})
