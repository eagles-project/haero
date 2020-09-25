include(UseLATEX)

add_latex_document(
  haero_design.tex
  INPUTS intro.tex
         physics.tex
         library.tex
         driver.tex
         api.tex
         driver_input.tex
         glossary.tex
  IMAGE_DIRS figures
  BIBFILES haero_design.bib
)
