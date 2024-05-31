#include "core/crd_transformations.hpp"
#include <cmath>

void dso::cartesian2spherical(double x, double y, double z, double &r,
                              double &glat, double &lon) noexcept {
  /* radius */
  r = std::sqrt(x * x + y * y + z * z);

  /* longitude */
  lon = std::atan2(y, x);

  /* geocentric latitude */
  glat = std::asin(z / r);

  return;
}
