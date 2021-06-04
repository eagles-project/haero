#include "toy_problem.hpp"
#include <sys/stat.h>
#include <sys/types.h>
#include <cstdio>

namespace haero {
namespace chemDriver {

void create_chem_files() {

  // Write out some test data to our current working directory.
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
  f = fopen("therm.dat", "w");
  fclose(f);

}

} // namespace chemDriver
} // namespace haero

