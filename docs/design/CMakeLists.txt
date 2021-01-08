include(UseLATEX)

# Generate a list of process sections and bibliographies.
file(GLOB PROCESS_SECTIONS
     LIST_DIRECTORIES FALSE
     RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
     processes/*.tex)

# Generate an appendix consisting of process-related sections. These sections
# have names that end with _app.tex and are stripped out of the previous list
# and collected into a single appendix.
set(PROCESS_APPENDIXES ${PROCESS_SECTIONS})
list(FILTER PROCESS_APPENDIXES INCLUDE REGEX .*_appendix.tex$)
set(PROC_APP_FILE ${CMAKE_CURRENT_BINARY_DIR}/process_appendix.tex)
file(WRITE ${PROC_APP_FILE} "\\chapter{Process-Specific Details}\n")
file(APPEND ${PROC_APP_FILE} "\\labelappendix{processes}\n\n")
file(APPEND ${PROC_APP_FILE}
     "This appendix contains information referenced in descriptions of various aerosol processes in \refchapter{processes}.\n\n")
foreach(app_file ${PROCESS_APPENDIXES})
  # Append the content of the appendix secton to our file.
  file(READ ${app_file} APP_TEXT)
  file(APPEND ${PROC_APP_FILE} ${APP_TEXT})

  # Remove this appendix from PROCESS_SECTIONS.
  list(REMOVE_ITEM PROCESS_SECTIONS ${app_file})
endforeach()

add_latex_document(
  haero_design.tex
  INPUTS intro.tex
         physics.tex
         library.tex
         processes.tex
         ${PROCESS_SECTIONS}
         testing.tex
         driver.tex
         api.tex
         driver_input.tex
         glossary.tex
         process_appendix.tex
  IMAGE_DIRS figures
  BIBFILES haero_design.bib
  TARGET_NAME design_doc
  EXCLUDE_FROM_ALL
)
set(DOC_TARGETS ${DOC_TARGETS} design_doc PARENT_SCOPE)
