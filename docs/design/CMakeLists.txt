include(UseLATEX)

add_latex_document(
  haero_design.tex
  INPUTS api.tex
         contexts.tex
         equations.tex
         input.tex
         intro.tex
         modes_and_species.tex
         tendencies.tex
  BIBFILES haero_design.bib
)
