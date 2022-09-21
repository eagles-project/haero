#include "kohler.hpp"

#include <iomanip>
#include <sstream>

namespace haero {
namespace kohler {

template <typename T>
std::string KohlerPolynomial<T>::mathematica_verification_program(const int n) const {
  std::ostringstream ss;
  /* mathematica doesn't like exponential notation, so we use setprecision
   * instead.*/
  ss << "kelvinCoeff = " << std::fixed << std::setprecision(16)
     << kelvin_a << ";\n";
  ss << "rhMin = " << rel_humidity_min << ";\n";
  ss << "rhMax = " << rel_humidity_max << ";\n";
  ss << "hygMin = " << hygro_min << ";\n";
  ss << "hygMax = " << hygro_max << ";\n";
  ss << "radMin = " << dry_radius_min_microns << ";\n";
  ss << "radMax = " << dry_radius_max_microns << ";\n";
  ss << "nn = " << n << ";\n";
  ss << "drh = (rhMax - rhMin)/(nn-1);\n";
  ss << "dhyg = (hygMax - hygMin)/(nn-1);\n";
  ss << "drad = (radMax - radMin)/(nn-1);\n";
  ss << "kinputs = Flatten[Table[{rhMin + i*drh, hygMin + j*dhyg, radMin + "
        "k*drad}, {i, 0, nn - 1}, {j, 0, nn - 1}, {k, 0, nn - 1}], 2];\n";
  ss << "kroots = Flatten[Table[NSolve[Log[kinputs[[i]][[1]]] r^4 + "
        "kelvinCoeff*(kinputs[[i]][[3]]^3 - r^3) + (kinputs[[i]][[2]] - "
        "Log[kinputs[[i]][[1]]]) kinputs[[i]][[3]]^3 r == 0 && r > 0, r, "
        "Reals], {i, Length[kinputs]}]];\n";
  ss << "kout = Table[r /. kroots[[i]], {i, Length[kroots]}];\n";
  ss << "(* Uncomment to overwrite data file: Export[\"" << HAERO_TEST_DATA_DIR
     << "/mm_kohler_roots.txt\", kout];*)\n";
  return ss.str();
}

template <typename T>
std::string KohlerPolynomial<T>::matlab_verification_program(const int n) const {
  std::ostringstream ss;
  ss << "clear; format long;\n";
  ss << "%% parameter bounds\n";
  ss << "kelvinCoeff = " << std::fixed << std::setprecision(16)
     << kelvin_a << ";\n";
  ss << "rhMin = " << rel_humidity_min << ";\n";
  ss << "rhMax = " << rel_humidity_max << ";\n";
  ss << "hygMin = " << hygro_min << ";\n";
  ss << "hygMax = " << hygro_max << ";\n";
  ss << "radMin = " << dry_radius_min_microns << ";\n";
  ss << "radMax = " << dry_radius_max_microns << ";\n";
  ss << "nn = " << n << ";\n";
  ss << "%% parameter inputs\n";
  ss << "relh = rhMin:(rhMax - rhMin)/(nn-1):rhMax;\n";
  ss << "hygro = hygMin:(hygMax-hygMin)/(nn-1):hygMax;\n";
  ss << "dry_rad = radMin:(radMax - radMin)/(nn - 1):radMax;\n";
  ss << "%% kohler polynomial companion matrix\n";
  ss << "cmat = zeros(4);\n";
  ss << "for i=1:3\n";
  ss << "  cmat(i+1,i) = 1;\n";
  ss << "end\n";
  ss << "%% solves: companion matrix eigenvalues\n";
  ss << "wet_rad_sol = zeros(1,nn^3);\n";
  ss << "idx = 1;\n";
  ss << "for i=1:nn\n";
  ss << "  logrh=log(relh(i));\n";
  ss << "  for j=1:nn\n";
  ss << "    hyg = hygro(j);\n";
  ss << "    for k=1:nn\n";
  ss << "      dradcube = dry_rad(k)^3;\n";
  ss << "      cmat(1,4) = -kelvinCoeff*dradcube/logrh;\n";
  ss << "      cmat(2,4) = (logrh-hyg)*dradcube/logrh;\n";
  ss << "      cmat(4,4) = kelvinCoeff/logrh;\n";
  ss << "      e = eig(cmat');\n";
  ss << "      wet_rad_sol(idx) = max(real(e));\n";
  ss << "      idx = idx + 1;\n";
  ss << "    end\n";
  ss << "  end\n";
  ss << "end\n";
  ss << "%% output: uncomment to overwrite data file\n";
  ss << "% writematrix(wet_rad_sol', \"" << HAERO_TEST_DATA_DIR
     << "/matlab_kohler_roots.txt\");\n";
  return ss.str();
}

// ETI
// double precison is required by the KohlerPolynomial class, so we intstantiate
// the two most common types here.
template struct KohlerPolynomial<ekat::Pack<double, HAERO_PACK_SIZE>>;
template struct KohlerPolynomial<double>;

}  // namespace kohler
}  // namespace haero
