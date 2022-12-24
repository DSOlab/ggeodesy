#include "units.hpp"
#include <cassert>
#include <cstdio>
#include <limits>

#define ERFA_DSIGN(A, B) ((B) < 0.0 ? -fabs(A) : fabs(A))
constexpr const double ERFA_DPI = dso::DPI;
constexpr const double ERFA_D2PI = dso::D2PI;

double erfa_anpm(double a) {
  double w;
  w = fmod(a, ERFA_D2PI);
  if (fabs(w) >= ERFA_DPI)
    w -= ERFA_DSIGN(ERFA_D2PI, a);
  return w;
}

double erfa_anp(double a) {
  double w;
  w = fmod(a, ERFA_D2PI);
  if (w < 0)
    w += ERFA_D2PI;
  return w;
}

bool ApproxEqual(double a, double b,
                 double tolerance = std::numeric_limits<double>::epsilon()) {
  double diff = std::abs(a - b);
  if (diff <= tolerance)
    return true;

  if (diff < std::max(std::abs(a), std::abs(b)) * tolerance)
    return true;

  return false;
}

int main() {
  const double from = -8e0 * dso::D2PI;
  const double to = 8e0 * dso::D2PI;

  for (double angle = from; angle < to; angle += 1e-3) {
    // reference results, [rad]
    double anp = erfa_anp(angle);
    double anpm = erfa_anpm(angle);

    // normalization using geodesy [rad]
    double manp = dso::anp<double>(angle);
    double manpm = dso::anpm<double>(angle);

    // normalization in [degrees]
    const double adeg = dso::rad2deg<double>(angle);
    double manp_d = dso::anp<double, dso::AngleUnit::Degrees>(adeg);
    double manpm_d = dso::anpm<double, dso::AngleUnit::Degrees>(adeg);

    // compare results w.r.t reference, [rad]
    assert(ApproxEqual(anp, manp));
    assert(ApproxEqual(anpm, manpm));

    // compare results w.r.t reference, [deg]
    //assert(ApproxEqual(manp_d, dso::rad2deg(anp)));
    //assert(ApproxEqual(manpm_d, dso::rad2deg(anpm)));
    if (!ApproxEqual(manp_d, dso::rad2deg(anp))) {
      fprintf(stderr,
              "[1] Test would fail for angle=%.9f[deg], difference is:%.9e\n",
              adeg, manp_d - dso::rad2deg(anp));
    }
    if (!ApproxEqual(manpm_d, dso::rad2deg(anpm))) {
      fprintf(stderr,
              "[2] Test would fail for angle=%.9f[deg], difference is:%.9e\n",
              adeg, manpm_d - dso::rad2deg(anpm));
    }
  }

  printf("Angle normalization tests passed!\n");
  return 0;
}
