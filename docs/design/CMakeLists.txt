include(UseLATEX)

# Generate a list of process sections and bibliographies.
file(GLOB PROCESS_SECTIONS
     LIST_DIRECTORIES FALSE
     RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
     processes/*.tex)
# TODO: Figure out per-section bibliographies
#file(GLOB PROCESS_BIBS
#     LIST_DIRECTORIES FALSE
#     RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
#     processes/*.bib)

add_latex_document(
  haero_design.tex
  INPUTS intro.tex
         physics.tex
         library.tex
         processes.tex
         ${PROCESS_SECTIONS}
         driver.tex
         api.tex
         driver_input.tex
         glossary.tex
  IMAGE_DIRS figures
  BIBFILES haero_design.bib #${PROCESS_BIBS}
)
