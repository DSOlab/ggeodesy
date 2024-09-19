#include "core/crd_transformations.hpp"

Eigen::Matrix<double, 3, 3> dso::geodetic2lvlh(double lat,
                                               double lon) noexcept {
  const double sphi = std::sin(lat);
  const double cphi = std::cos(lat);
  const double slmb = std::sin(lon);
  const double clmb = std::cos(lon);

  Eigen::Matrix<double, 3, 3> enu;
  /* unit vector e, i.e. east */
  enu(0, 0) = -slmb;
  enu(1, 0) = clmb;
  enu(2, 0) = 0e0;
  /* unit vector n, i.e. north */
  enu(0, 1) = -sphi * clmb;
  enu(1, 1) = -sphi * slmb;
  enu(2, 1) = cphi;
  /* unit vector u, i.e. up */
  enu(0, 2) = cphi * clmb;
  enu(1, 2) = cphi * slmb;
  enu(2, 2) = sphi;

  return enu;
}
