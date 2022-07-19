#include <cstdio>

#include "catch2/catch.hpp"
#include "haero/processes/vehkamaki2002.hpp"

using namespace haero;

// These tests aren't really unit tests in the usual sense--they are "smoke
// tests" in that they exercise the parameterizations of Vehkamaki et al
// (2002) and generate Python scripts that generate plots. These plots can be
// used to compare the results of the parameterizations with the calculations
// from classical nucleation theory (e.g. given in Chapter 11 of the Third
// Edition of Seinfeld and Pandis, 2016).

namespace {

// Create a vector of values between two endpoints v1 and v2 with a step size
// dv.
std::vector<PackType> create_inputs(Real v1, Real v2, Real dv) {
  size_t len = 1 + static_cast<size_t>(std::ceil((v2 - v1) / dv));
  std::vector<PackType> v(len);
  for (int i = 0; i < len; ++i) {
    v[i] = v1 + i * dv;
  }
  return v;
}

// Create a vector of output generated from traversing an input array and
// applying the given function to it.
std::vector<PackType> generate_output(const std::vector<PackType>& x,
                                      std::function<PackType(PackType)> f) {
  std::vector<PackType> y(x.size());
  for (int i = 0; i < x.size(); ++i) {
    y[i] = f(x[i]);
  }
  return y;
}

// Writes a Python script that generates a plot of a curve y(x).
void generate_curve_plot(const std::string& filename, const std::string& title,
                         const std::string& x_label, const std::string& x_scale,
                         const std::string& y_label, const std::string& y_scale,
                         const std::vector<PackType>& x,
                         const std::vector<PackType>& y, Real xmin = 1.0,
                         Real xmax = -1.0, Real ymin = 1.0, Real ymax = -1.0) {
  FILE* fp = fopen(filename.c_str(), "w");
  fprintf(fp, "import matplotlib.pyplot as plt\n\n");

  fprintf(fp, "# Plot data\n");
  fprintf(fp, "x = [");
  for (int i = 0; i < x.size(); ++i) {
    fprintf(fp, "%g, ", x[i][0]);
  }
  fprintf(fp, "]\n");
  fprintf(fp, "y = [");
  for (int i = 0; i < y.size(); ++i) {
    fprintf(fp, "%g, ", y[i][0]);
  }
  fprintf(fp, "]\n\n");

  fprintf(fp, "fig, ax = plt.subplots()\n");
  fprintf(fp, "ax.plot(x, y)\n\n");
  fprintf(fp,
          "ax.set(xlabel='%s', xscale='%s', ylabel='%s', yscale='%s', "
          "title='%s')\n",
          x_label.c_str(), x_scale.c_str(), y_label.c_str(), y_scale.c_str(),
          title.c_str());
  if (xmin < xmax) {  // if x limits are specified
    fprintf(fp, "ax.set_xlim(%g, %g)\n", xmin, xmax);
  }
  if (ymin < ymax) {  // if y limits are specified
    fprintf(fp, "ax.set_ylim(%g, %g)\n", ymin, ymax);
  }
  fprintf(fp, "ax.grid()\n\n");

  size_t dot = filename.find(".");
  std::string plotfile = filename.substr(0, dot) + ".png";
  fprintf(fp, "print('Generating plot: %s')\n", plotfile.c_str());
  fprintf(fp, "fig.savefig('%s')\n", plotfile.c_str());
  fclose(fp);
}

// Create a vector of output generated from traversing two input arrays and
// applying the given function to it.
std::vector<std::vector<PackType>> generate_output(
    const std::vector<PackType>& x, const std::vector<PackType>& y,
    std::function<PackType(PackType, PackType)> f) {
  std::vector<std::vector<PackType>> z(x.size());
  for (int i = 0; i < x.size(); ++i) {
    z[i].resize(y.size());
    for (int j = 0; j < y.size(); ++j) {
      z[i][j] = f(x[i], y[j]);
    }
  }
  return z;
}

// Writes a Python script that generates a plot of a family of curves y(x)
// parameterized by the values in p.
void generate_multicurve_plot(const std::string& filename,
                              const std::string& title,
                              const std::string& x_label,
                              const std::string& y_label,
                              const std::string& p_label,
                              const std::vector<PackType>& x,
                              const std::vector<std::vector<PackType>>& y,
                              const std::vector<PackType>& p) {
  FILE* fp = fopen(filename.c_str(), "w");
  fprintf(fp, "import matplotlib.pyplot as plt\n\n");

  fprintf(fp, "# Plot data\n");
  fprintf(fp, "x = [");
  for (int i = 0; i < x.size(); ++i) {
    fprintf(fp, "%g, ", x[i][0]);
  }
  fprintf(fp, "]\n");
  for (int i = 0; i < y.size(); ++i) {
    fprintf(fp, "y%d = [", i + 1);
    for (int j = 0; j < y[i].size(); ++j) {
      fprintf(fp, "%g, ", y[i][j][0]);
    }
    fprintf(fp, "]\n");
  }
  fprintf(fp, "p = [");
  for (int i = 0; i < p.size(); ++i) {
    fprintf(fp, "%g, ", p[i][0]);
  }
  fprintf(fp, "]\n\n");

  for (int i = 0; i < y.size(); ++i) {
    fprintf(fp, "plt.plot(x, y%d, label='%s = %g')\n", i + 1, p_label.c_str(),
            p[i][0]);
  }
  fprintf(fp, "plt.legend()\n");
  fprintf(fp, "plt.title('%s')\n", title.c_str());
  fprintf(fp, "plt.xlabel('%s')\n", x_label.c_str());
  fprintf(fp, "plt.ylabel('%s')\n\n", y_label.c_str());

  size_t dot = filename.find(".");
  std::string plotfile = filename.substr(0, dot) + ".png";
  fprintf(fp, "print('Generating plot: %s')\n", plotfile.c_str());
  fprintf(fp, "plt.savefig('%s')\n", plotfile.c_str());
  fclose(fp);
}

}  // anonymous namespace

// Compare the output of this test with Vehkamaki et al (2002), figure 7.
TEST_CASE("vehkamaki2002_figure_7") {
  auto temp = create_inputs(190.15, 290.15, 5);
  auto rel_hum = create_inputs(0.05, 1, 0.05);
  auto J = generate_output(temp, rel_hum, [=](PackType T, PackType RH) {
    return vehkamaki2002::h2so4_nucleation_threshold(T, RH);
  });
  generate_multicurve_plot("vehkamaki2002_figure_7.py", "Nucleation rate",
                           "Relative humidity", "Nucleation rate [#/cm3]",
                           "T [K]", rel_hum, J, temp);
}

// Compare the output of this test with Seinfeld and Pandis, figure 11.11.
// Alas, we can only test one temperature value within the region of validity!
TEST_CASE("vehkamaki2002_sp_figure_11.11") {
  PackType temp(273.0);
  auto rel_hum = create_inputs(0.2, 1.0, 0.1);
  auto x_crit = generate_output(rel_hum, [=](PackType RH) {
    auto c_h2so4 = vehkamaki2002::h2so4_nucleation_threshold(temp, RH);
    return vehkamaki2002::h2so4_critical_mole_fraction(c_h2so4, temp, RH);
  });
  generate_curve_plot("vehkamaki2002_sp_figure_11.11.py",
                      "Critical cluster mole fraction", "Relative humidity",
                      "linear", "Mole fraction", "linear", rel_hum, x_crit);
}

// Compare the output of this test with Vehkamaki et al (2002), figure 8
// (dotted line).
TEST_CASE("vehkamaki2002_figure_8_dotted") {
  PackType temp(236.0), rel_hum(0.55);
  auto c_h2so4 = create_inputs(1e6, 100e7, 1e7);
  auto J = generate_output(c_h2so4, [=](PackType c) {
    auto x_crit = vehkamaki2002::h2so4_critical_mole_fraction(c, temp, rel_hum);
    return vehkamaki2002::nucleation_rate(c, temp, rel_hum, x_crit);
  });
  generate_curve_plot("vehkamaki2002_figure_8_dotted.py", "Nucleation Rate",
                      "H2SO4 concentration", "log", "J [#/cm3]", "log", c_h2so4,
                      J);
}

// Compare the output of this test with Vehkamaki et al (2002), figure 8
// (dashed line).
TEST_CASE("vehkamaki2002_figure_8_dashed") {
  PackType temp(236.0), rel_hum(0.55);
  auto c_h2so4 = create_inputs(1e6, 100e7, 1e7);
  auto n_tot = generate_output(c_h2so4, [=](PackType c) {
    auto x_crit = vehkamaki2002::h2so4_critical_mole_fraction(c, temp, rel_hum);
    return vehkamaki2002::num_critical_molecules(c, temp, rel_hum, x_crit);
  });
  generate_curve_plot("vehkamaki2002_figure_8_dashed.py",
                      "Critical cluster size", "H2SO4 concentration", "log",
                      "n_tot", "linear", c_h2so4, n_tot);
}

// Compare the output of this test with Vehkamaki et al (2002), figure 9.
TEST_CASE("vehkamaki2002_figure_9") {
  PackType temp(298.0), rel_hum = 0.382;
  auto c_h2so4 = create_inputs(7.5e9, 3.5e10, 1e8);
  auto J = generate_output(c_h2so4, [=](PackType c) {
    auto x_crit = vehkamaki2002::h2so4_critical_mole_fraction(c, temp, rel_hum);
    return vehkamaki2002::nucleation_rate(c, temp, rel_hum, x_crit);
  });
  generate_curve_plot("vehkamaki2002_figure_9a.py", "298K, RH=38.2%",
                      "H2SO4 concentration [#/cc]", "linear",
                      "Nucleation rate [#/cc/s]", "log", c_h2so4, J, 0, 4e10,
                      1e-3, 1e6);
  rel_hum = 0.523;
  J = generate_output(c_h2so4, [=](PackType c) {
    auto x_crit = vehkamaki2002::h2so4_critical_mole_fraction(c, temp, rel_hum);
    return vehkamaki2002::nucleation_rate(c, temp, rel_hum, x_crit);
  });
  generate_curve_plot("vehkamaki2002_figure_9b.py", "298K, RH=52.3%",
                      "H2SO4 concentration [#/cc]", "linear",
                      "Nucleation rate [#/cc/s]", "log", c_h2so4, J, 0, 4e10,
                      1e-3, 1e6);
}

// Compare the output of this test with Vehkamaki et al (2002), figure 10.
TEST_CASE("vehkamaki2002_figure_10") {
  PackType temp(295.15), rel_hum = 0.075;
  auto c_h2so4 = create_inputs(1e10, 1e12, 1e9);
  auto J = generate_output(c_h2so4, [=](PackType c) {
    auto x_crit = vehkamaki2002::h2so4_critical_mole_fraction(c, temp, rel_hum);
    return vehkamaki2002::nucleation_rate(c, temp, rel_hum, x_crit);
  });
  generate_curve_plot("vehkamaki2002_figure_10a.py", "295.15K, RH=7.5%",
                      "H2SO4 concentration [#/cc]", "log",
                      "Nucleation rate [#/cc/s]", "log", c_h2so4, J, 1e10, 1e12,
                      1e-2, 1e5);
  rel_hum = 0.153;
  J = generate_output(c_h2so4, [=](PackType c) {
    auto x_crit = vehkamaki2002::h2so4_critical_mole_fraction(c, temp, rel_hum);
    return vehkamaki2002::nucleation_rate(c, temp, rel_hum, x_crit);
  });
  generate_curve_plot("vehkamaki2002_figure_10b.py", "295.15K, RH=15.3%",
                      "H2SO4 concentration [#/cc]", "log",
                      "Nucleation rate [#/cc/s]", "log", c_h2so4, J, 1e10, 1e12,
                      1e-2, 1e5);
}
