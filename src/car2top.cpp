#include "geodesy.hpp"
#include <cstdio>

void dso::core::dcar2top(double xi, double yi, double zi, double dx, double dy,
                         double dz, double semi_major, double flattening,
                         double &east, double &north, double &up) noexcept {

  // Ellipsoidal coordinates of reference point.
  double phi_i, lambda_i, h_i;

  // Cartesian to ellipsoidal for reference point.
  core::car2ell(xi, yi, zi, semi_major, flattening, phi_i, lambda_i, h_i);

  // Trigonometric numbers.
  const double cosf{std::cos(phi_i)};
  const double cosl{std::cos(lambda_i)};
  const double sinf{std::sin(phi_i)};
  const double sinl{std::sin(lambda_i)};

  // Topocentric vector.
  north = -sinf * cosl * dx - sinf * sinl * dy + cosf * dz;
  east = -sinl * dx + cosl * dy;
  up = cosf * cosl * dx + cosf * sinl * dy + sinf * dz;

  // Finished.
  return;
}

Eigen::Matrix<double, 3, 1>
dso::core::dcar2top(const Eigen::Matrix<double, 3, 1> &r,
                    const Eigen::Matrix<double, 3, 1> &dr, double semi_major,
                    double flattening) noexcept {

  // Cartesian to ellipsoidal for reference point.
  const Eigen::Matrix<double, 3, 1> lfh =
      dso::core::car2ell(r, semi_major, flattening);
  return dso::topocentric_matrix(lfh) * dr;
}
