/// @file ell2car.cpp
/// @brief Transform ellipsoidal to cartesian coordinates.
#include "geodesy.hpp"

/// @brief Ellipsoidal to cartesian coordinates.
///
/// Transform (geocentric) cartesian coordinates (on the ellipsoid) to
/// ellipsoidal coordinates. Units are meters and radians.
///
/// @param[in]   phi    Ellipsoidal latitude (radians)
/// @param[in]   lambda Ellipsoidal longtitude (radians)
/// @param[in]   h      Ellipsoidal height (meters)
/// @param[in]   e      The reference ellipsoid (dso::Ellipsoid)
/// @param[out]  x      Cartesian x-component (meters)
/// @param[out]  y      Cartesian y-component (meters)
/// @param[out]  z      Cartesian z-component (meters)
/// @throw              Does not throw.
void ell2car(double lambda, double phi, double h, const dso::Ellipsoid &e, double &x,
             double &y, double &z) noexcept {
  // Eccentricity squared.
  double e2{e.eccentricity_squared()};

  // Radius of curvature in the prime vertical.
  const double N{e.N(phi)};

  // Trigonometric numbers.
  const double sinf{std::sin(phi)};
  const double cosf{std::cos(phi)};
  const double sinl{std::sin(lambda)};
  const double cosl{std::cos(lambda)};

  // Compute geocentric rectangular coordinates.
  x = (N + h) * cosf * cosl;
  y = (N + h) * cosf * sinl;
  z = ((1e0 - e2) * N + h) * sinf;

  // Finished.
  return;
}

Eigen::Matrix<double,3,1> dso::ell2car(const Eigen::Matrix<double,3,1> &lfh, const dso::Ellipsoid &e) noexcept {
  // Eccentricity squared.
  double e2{e.eccentricity_squared()};

  // Radius of curvature in the prime vertical.
  const double N{e.N(lfh(1))};
  
  // Trigonometric numbers.
  const double sinf{std::sin(lfh(1))};
  const double cosf{std::cos(lfh(1))};
  const double sinl{std::sin(lfh(0))};
  const double cosl{std::cos(lfh(0))};

  // Compute geocentric rectangular coordinates.
  const double x = (N + lfh(2)) * cosf * cosl;
  const double y = (N + lfh(2)) * cosf * sinl;
  const double z = ((1e0 - e2) * N + lfh(2)) * sinf;

  // Finished.
  /*const*/ double data[] = {x,y,z};
  return Eigen::Map<Eigen::Matrix<double,3,1>>(data,3);
}
