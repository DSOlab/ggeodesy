/// @file ellipsoid_core.hpp
/// @brief Core function of ellipsoidal geometry
/// [1] H. Moritz, GEODETIC REFERENCE SYSTEM 1980,
/// https://geodesy.geology.ohio-state.edu/course/refpapers/00740128.pdf

#ifndef __ELLIPSOID_GEOMETRY_CORE_HPP__
#define __ELLIPSOID_GEOMETRY_CORE_HPP__

#include <cmath>

namespace dso {

/// @brief core namespace holds the core of ellipsoid-related functions
namespace core {

/// @brief Compute the squared eccentricity.
///
/// Compute the squared (first) eccentricity (i.e. \f$e^2\f$) given the
/// flattening of an ellipsoid, aka \f$ e^2 = \frac{a^2-b^2}{a^2} = (2-f)*f \f$
/// @param[in] f flattening
/// @return      squared eccentricity
inline constexpr double eccentricity_squared(double f) noexcept {
  return (2e0 - f) * f;
}

/// @brief Compute the third flattening
///
/// Compute the third flattening (usually denoted as \f$ n \f$), given the
/// flattening \f$ f \f$, aka \f$ n = \frac{a-b}{a+b} = f/(2-f) \f$.
/// @param[in] f flattening
/// @return third flattening \f$ n \f$
/// Reference [2]
inline constexpr double third_flattening(double f) noexcept {
  return f / (2e0 - f);
}

/// @brief Compute the semi-minor axis (aka \f$b\f$)
///
/// Compute the semi-minor axis of an ellipsoid(i.e. \f$e^2\f$), given the
/// flattening and the semi-major axis, aka \f$ \beta = \alpha * (1-f) \f$
/// @param[in] f flattening
/// @param[in] a semi-major axis (meters)
/// @return      semi-minor axis (meters)
inline constexpr double semi_minor(double a, double f) noexcept {
  return a * (1e0 - f);
}

/// @brief the linear eccentricity
///
/// Compute the linear eccentricity of an ellipsoid (i.e. E) from the formula
/// \f$ E = \sqrt{a^2 - b^2} \f$
/// @param[in] f flattening
/// @param[in] a semi-major axis (meters)
/// @return      semi-minor axis (meters)
inline
#if defined(__GNUC__) && !defined(__llvm__)
    constexpr
#endif
    double
    linear_eccentricity(double a, double f) noexcept {
  const double b = semi_minor(a, f);
  return std::sqrt(a * a - b * b);
}

/// @brief Polar radius of curvature
///
/// Compute the polar radius of curvature of an ellipsoid (i.e. c) from the
/// formula \f$ c = a^2 / b\f$
/// @param[in] f flattening
/// @param[in] a semi-major axis (meters)
/// @return      semi-minor axis (meters)
inline constexpr double polar_radius_of_curvature(double a, double f) noexcept {
  const double b = semi_minor(a, f);
  return a * a / b;
}

/// @brief Compute the normal radius of curvature at a given latitude (on
///        a reference ellipsoid).
///
/// @param[in] lat The latitude (radians)
/// @param[in] a   The ellipsoid's semi-major axis (meters)
/// @param[in] b   The ellipsoid's semi-minor axis (meters)
/// @return        The normal radius of curvature (meters)
/// @throw         Does not throw (see notes)
///
/// @note If the denominator (den) is zero then funny things could happen;
///       this however should **never** occur for any reference ellipsoid.
///
/// @see "Physical Geodesy", pg. 194
/// @see https://en.wikipedia.org/wiki/Earth_radius
inline
#if defined(__GNUC__) && !defined(__llvm__)
    constexpr
#endif
    double
    N(double lat, double a, double b) noexcept {
  const double cosf{std::cos(lat)};
  const double sinf{std::sin(lat)};
  const double acosf{a * cosf};
  const double bsinf{b * sinf};
  const double den{std::sqrt(acosf * acosf + bsinf * bsinf)};
  return (a * a) / den;
}

/// @brief Compute the meridional radii of curvature at a given latitude
///        on a reference ellipsoid.
///
/// @param[in] lat The latitude in (radians)
/// @param[in] a   The ellipsoid's semi-major axis (meters)
/// @param[in] b   The ellipsoid's semi-minor axis (meters)
/// @return        The meridional radius of curvature (meters)
/// @throw         Does not throw (see notes)
///
/// @see https://en.wikipedia.org/wiki/Earth_radius
///
inline
#if defined(__GNUC__) && !defined(__llvm__)
    constexpr
#endif
    double
    M(double lat, double a, double b) noexcept {
  const double cosf{std::cos(lat)};
  const double sinf{std::sin(lat)};
  const double acosf{a * cosf};
  const double bsinf{b * sinf};
  const double tmpd{acosf * acosf + bsinf * bsinf};
  return ((a * b) / tmpd) * ((a * b) / std::sqrt(tmpd));
}

/// @brief Compute the geocentric latitude
///
/// The geocentric latitude is the angle between the equatorial plane and the
/// radius from the centre to a point on the surface. The relation between the
/// geocentric latitude(\f$\theta\f$) and the geodetic latitude(\f$\phi\f$) is
/// \f$ \theta (\phi) = tan^{-1} ((1-f)^2 tan(\phi)) \f$
/// The geodetic and geocentric latitudes are equal at the equator and at the
/// poles but at other latitudes they differ by a few minutes of arc.
/// 
/// @param[in] lat The (geodetic) latitude in [rad]
/// @param[in] f   The ellipsoid's flattening [-]
/// @return        The geocentric latitude [rad]
///
/// Reference Torge, 2001, Eq. 4.11
inline
#if defined(__GNUC__) && !defined(__llvm__)
    constexpr
#endif
    double
    geocentric_latitude(double lat, double f) noexcept {
  return std::atan((1e0 - f) * (1e0 - f) * std::tan(lat));
}

/// @brief Compute the parametric or reduced latitude
///
/// The parametric or reduced latitude, \f$ beta \f$ is defined by the radius
/// drawn from the centre of the ellipsoid to that point Q on the surrounding
/// sphere (of radius a) which is the projection parallel to the Earth's axis
/// of a point P on the ellipsoid at latitude \f$ \phi \f$
/// @param[in] lat The (geodetic) latitude in radians
/// @param[in] f   The ellipsoid's flattening
/// @return        The parametric or reduced latitude at lat in radians
///
/// Reference Torge, 2001, Eq. 4.11
inline
#if defined(__GNUC__) && !defined(__llvm__)
    constexpr
#endif
    double
    reduced_latitude(double lat, double f) noexcept {
  return std::atan((1e0 - f) * std::tan(lat));
}

} // core

} // dso

#endif
