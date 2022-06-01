/// @file car2ell.cpp
/// @brief Transformation of ellipsoidal to cartesian coordinates.

#include "geodesy.hpp"

using dso::VECTOR3;

/// @brief Cartesian to ellipsoidal.
///
/// Transform cartesian, geocentric coordinates (x, y, z) to ellipsoidal (i.e.
/// latitude, longtitude, ellispoidal height). All units are meters and
/// radians.
///
/// @param[in]  x      Cartesian x-component (meters)
/// @param[in]  y      Cartesian y-component (meters)
/// @param[in]  z      Cartesian z-component (meters)
/// @param[in]  semi_major   Semi-major axis of the (reference) ellipsoid
//                     (meters)
/// @param[in]  flattening   Flattening of the (reference) ellipsoid
/// @param[out] phi    Ellipsoidal latitude (radians)
/// @param[out] lambda Ellipsoidal longtitude (radians)
/// @param[out] h      Ellipsoidal height (meters)
/// @throw             Does not throw.
///
/// @see Fukushima, T., "Transformation from Cartesian to geodetic coordinates
///      accelerated by Halley's method", J. Geodesy (2006), 79(12): 689-693
///
void dso::core::car2ell(double x, double y, double z, double semi_major, double flattening,
             double &phi, double &lambda, double &h) noexcept {
  // Functions of ellipsoid parameters.
  const double aeps2{semi_major * semi_major * 1e-32};
  const double e2{(2.0e0 - flattening) * flattening};
  const double e4t{e2 * e2 * 1.5e0};
  const double ep2{1.0e0 - e2};
  const double ep{std::sqrt(ep2)};
  const double aep{semi_major * ep};

  // Compute Coefficients of (Modified) Quartic Equation
  // Remark: Coefficients are rescaled by dividing by 'a'

  // Compute distance from polar axis squared.
  double p2{x * x + y * y};

  // Compute longitude lambda.
  if (p2) {
    lambda = std::atan2(y, x);
  } else {
    lambda = .0e0;
  }

  // Ensure that Z-coordinate is unsigned.
  double absz{std::abs(z)};

  if (p2 > aeps2) { // Continue unless at the poles
    // Compute distance from polar axis.
    double p{std::sqrt(p2)};
    // Normalize.
    double s0{absz / semi_major};
    double pn{p / semi_major};
    double zp{ep * s0};
    // Prepare Newton correction factors.
    double c0{ep * pn};
    double c02{c0 * c0};
    double c03{c02 * c0};
    double s02{s0 * s0};
    double s03{s02 * s0};
    double a02{c02 + s02};
    double a0{std::sqrt(a02)};
    double a03{a02 * a0};
    double d0{zp * a03 + e2 * s03};
    double f0{pn * a03 - e2 * c03};
    // Prepare Halley correction factor.
    double b0{e4t * s02 * c02 * pn * (a0 - ep)};
    double s1{d0 * f0 - b0 * s0};
    double cp{ep * (f0 * f0 - b0 * c0)};
    // Evaluate latitude and height.
    phi = ::atan(s1 / cp);
    double s12{s1 * s1};
    double cp2{cp * cp};
    h = (p * cp + absz * s1 - semi_major * std::sqrt(ep2 * s12 + cp2)) /
        std::sqrt(s12 + cp2);
  } else { // Special case: pole.
    phi = dso::DPI / 2e0;
    h = absz - aep;
  }

  // Restore sign of latitude.
  if (z < 0.e0)
    phi = -phi;

  // Finished.
  return;
}

VECTOR3 dso::car2ell(const VECTOR3 &xyz, double semi_major,
                     double flattening) noexcept {
  const double x = xyz(0);
  const double y = xyz(1);
  const double z = xyz(2);

  // Functions of ellipsoid parameters.
  const double aeps2{semi_major * semi_major * 1e-32};
  const double e2{(2.0e0 - flattening) * flattening};
  const double e4t{e2 * e2 * 1.5e0};
  const double ep2{1.0e0 - e2};
  const double ep{std::sqrt(ep2)};
  const double aep{semi_major * ep};

  // Compute Coefficients of (Modified) Quartic Equation
  // Remark: Coefficients are rescaled by dividing by 'a'

  // Compute distance from polar axis squared.
  double p2{x * x + y * y};

  double lambda, phi, hgt;

  // Compute longitude lambda.
  if (p2) {
    lambda = std::atan2(y, x);
  } else {
    lambda = .0e0;
  }

  // Ensure that Z-coordinate is unsigned.
  double absz{std::abs(z)};

  if (p2 > aeps2) { // Continue unless at the poles
    // Compute distance from polar axis.
    const double p{std::sqrt(p2)};
    // Normalize.
    const double s0{absz / semi_major};
    const double pn{p / semi_major};
    const double zp{ep * s0};
    // Prepare Newton correction factors.
    const double c0{ep * pn};
    const double c02{c0 * c0};
    const double c03{c02 * c0};
    const double s02{s0 * s0};
    const double s03{s02 * s0};
    const double a02{c02 + s02};
    const double a0{std::sqrt(a02)};
    const double a03{a02 * a0};
    const double d0{zp * a03 + e2 * s03};
    const double f0{pn * a03 - e2 * c03};
    // Prepare Halley correction factor.
    const double b0{e4t * s02 * c02 * pn * (a0 - ep)};
    const double s1{d0 * f0 - b0 * s0};
    const double cp{ep * (f0 * f0 - b0 * c0)};
    // Evaluate latitude and height.
    phi = ::atan(s1 / cp);
    const double s12{s1 * s1};
    const double cp2{cp * cp};
    hgt = (p * cp + absz * s1 - semi_major * std::sqrt(ep2 * s12 + cp2)) /
          std::sqrt(s12 + cp2);
  } else { // Special case: pole.
    phi = dso::DPI / 2e0;
    hgt = absz - aep;
  }

  // Restore sign of latitude.
  if (z < 0.e0)
    phi = -phi;

  // Finished.
  #ifdef USE_EIGEN
  double data[] = {lambda, phi, hgt};
  return Eigen::Map<VECTOR3>(data,3);
#else
  return Vector3({lambda, phi, hgt});
#endif
}

/// @brief Cartesian to ellipsoidal.
///
/// @tparam     E      The reference ellipsoid (i.e. one of dso::ellipsoid).
/// @param[in]  x      Cartesian x-component, meters.
/// @param[in]  y      Cartesian y-component, meters.
/// @param[in]  z      Cartesian z-component, meters.
/// @param[out] phi    Ellipsoidal latitude, radians.
/// @param[out] lambda Ellipsoidal longtitude, radians.
/// @param[out] h      Ellipsoidal height, meters.
///
/// @see  dso::core::car2ell
//template <ellipsoid E>
//void car2ell(double x, double y, double z, double &phi, double &lambda,
//             double &h) noexcept {
//  constexpr double semi_major{ellipsoid_traits<E>::a};
//  constexpr double flattening{ellipsoid_traits<E>::f};
//  core::car2ell(x, y, z, semi_major, flattening, phi, lambda, h);
//  return;
//}

/// @brief Cartesian to ellipsoidal.
///
/// @param[in]  x      Cartesian x-component, meters.
/// @param[in]  y      Cartesian y-component, meters.
/// @param[in]  z      Cartesian z-component, meters.
/// @param[in]  e      An dso::Ellipsoid instance
/// @param[out] phi    Ellipsoidal latitude, radians.
/// @param[out] lambda Ellipsoidal longtitude, radians.
/// @param[out] h      Ellipsoidal height, meters.
///
/// @see  dso::core::car2ell
//void car2ell(double x, double y, double z, const Ellipsoid &e, double &phi,
//             double &lambda, double &h) noexcept {
//  double semi_major{e.semi_major()};
//  double flattening{e.flattening()};
//  core::car2ell(x, y, z, semi_major, flattening, phi, lambda, h);
//  return;
//}

/// @brief Cartesian to ellipsoidal.
///
/// @param[in]  x      Cartesian x-component, meters.
/// @param[in]  y      Cartesian y-component, meters.
/// @param[in]  z      Cartesian z-component, meters.
/// @param[in]  e      The reference ellipsoid (i.e. one of dso::ellipsoid).
/// @param[out] phi    Ellipsoidal latitude, radians.
/// @param[out] lambda Ellipsoidal longtitude, radians.
/// @param[out] h      Ellipsoidal height, meters.
///
/// @see  dso::core::car2ell
//void car2ell(double x, double y, double z, ellipsoid e, double &phi,
//             double &lambda, double &h) noexcept {
//  car2ell(x, y, z, Ellipsoid(e), phi, lambda, h);
//}
