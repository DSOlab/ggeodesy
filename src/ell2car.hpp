/// @file ell2car.hpp
/// @brief Transform ellipsoidal to cartesian coordinates.

#ifndef __ELLIPSOIDAL_TO_CARTESIAN__
#define __ELLIPSOIDAL_TO_CARTESIAN__

#include "ellipsoid.hpp"
#include <cmath>

namespace dso {

/// @brief Ellipsoidal to cartesian coordinates.
///
/// Transform (geocentric) cartesian coordinates (on the ellipsoid) to
/// ellipsoidal coordinates. Units are meters and radians.
///
/// @tparam      E      The reference ellipsoid (i.e. one of dso::ellipsoid).
/// @param[in]   phi    Ellipsoidal latitude (radians)
/// @param[in]   lambda Ellipsoidal longtitude (radians)
/// @param[in]   h      Ellipsoidal height (meters)
/// @param[out]  x      Cartesian x-component (meters)
/// @param[out]  y      Cartesian y-component (meters)
/// @param[out]  z      Cartesian z-component (meters)
/// @throw              Does not throw.
///
template <ellipsoid E>
void ell2car(double phi, double lambda, double h, double &x, double &y,
             double &z) noexcept {
  // Eccentricity squared.
  constexpr double e2{dso::eccentricity_squared<E>()};

  // Radius of curvature in the prime vertical.
  const double N{dso::N<E>(phi)};

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
void ell2car(double phi, double lambda, double h, const Ellipsoid &e, double &x,
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
  z = ((1.0e0 - e2) * N + h) * sinf;

  // Finished.
  return;
}

} // namespace dso

#endif
