#include "geodesy.hpp"

using dso::MATRIX3x3;
using dso::VECTOR3;

void dso::core::dcar2top(double xi, double yi, double zi, double dx, double dy, double dz,
              double semi_major, double flattening, double &east, double &north,
              double &up) noexcept {

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

VECTOR3 dso::core::dcar2top(const VECTOR3 &r, const VECTOR3 &dr, double semi_major,
                 double flattening) noexcept {

  // Cartesian to ellipsoidal for reference point.
  const VECTOR3 lfh = dso::car2ell(r, semi_major, flattening);

  // Trigonometric numbers.
  //const double cosf{std::cos(lfh.y())};
  //const double sinf{std::sin(lfh.y())};
  //const double sinl{std::sin(lfh.x())};
  //const double cosl{std::cos(lfh.x())};

  //// Topocentric vector.
  //const double north =
  //    -sinf * cosl * dr.x() - sinf * sinl * dr.y() + cosf * dr.z();
  //const double east = -sinl * dr.x() + cosl * dr.y();
  //const double up = cosf * cosl * dr.x() + cosf * sinl * dr.y() + sinf * dr.z();

  // Finished.
  // return Vector3({east, north, up});

  return dso::topocentric_matrix(lfh) * dr;
}
