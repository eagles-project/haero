#include "toy_problem.hpp"

#include <sys/stat.h>
#include <sys/types.h>

#include <cstdio>

namespace haero {
namespace chem_driver {

void create_chem_files() {
  // create the TChem input file required for the toy problem
  // Note: the goal is to ultimately write this given the information from the
  // input yaml file, but I need to investigate whether TChem has moved to a
  // yaml input spec that might make this easier
  const char* chem_inp = R"INPUT(ELEMENTS
X /1/
END
SPECIES
X X2
END
THERM ALL
    300.000  1000.000  5000.000
X                        X  1               G   200.000  6000.000 1000.00      1
 2.50000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00    2
-7.45375000E+02 4.37967000E+00 2.50000000E+00 0.00000000E+00 0.00000000E+00    3
 0.00000000E+00 0.00000000E+00-7.45375000E+02 4.37967000E+00                   4
X2                       X  2               G   200.000  6000.000 1000.00      1
 2.50000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00    2
-7.45375000E+02 4.37967000E+00 2.50000000E+00 0.00000000E+00 0.00000000E+00    3
 0.00000000E+00 0.00000000E+00-7.45375000E+02 4.37967000E+00                   4
END
REACTIONS
X2=>2X      1E+0    1   1
X+X=>X2      1E+0    1   1
END
  )INPUT";
  FILE* f = fopen("chem.inp", "w");
  fprintf(f, "%s", chem_inp);
  fclose(f);
  // Note: as things are currently working, reaction rates are passed via the
  // input yaml file, whereas TChem initially looks in therm.dat.
  // As such, we need to create this empty file
  f = fopen("therm.dat", "w");
  fclose(f);
}

}  // namespace chem_driver
}  // namespace haero
