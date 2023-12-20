#include "crd_transformations.hpp"
#include <cmath>

void dso::spherical2cartesian(double r, double glat, double lon, double &x,
                         double &y, double &z) noexcept {

  const double sf = std::sin(glat);
  const double cf = std::cos(glat);
  const double sl = std::sin(lon);
  const double cl = std::cos(lon);

  x = r * cf * cl;
  y = r * cf * sl;
  z = r * sf;

  return;
}
