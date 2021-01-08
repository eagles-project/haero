
# Design document.
find_package(LATEX COMPONENTS PDFLATEX BIBTEX)
if (LATEX_PDFLATEX_FOUND AND LATEX_BIBTEX_FOUND)
  add_subdirectory(design)
endif()

# Doxygen reference documentation
find_package(Doxygen OPTIONAL_COMPONENTS dot QUIET)
if (DOXYGEN_FOUND)
  add_subdirectory(doxygen)
endif()

# Create a "docs" target that builds all documentation.
add_custom_target(docs DEPENDS ${DOC_TARGETS})
