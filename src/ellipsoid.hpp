///
/// @file ellipsoid.hpp
///
/// @brief Definition and basic setup of an ellipsoid (kinda) class.
///
/// This file defines a list of frequently used (referennce) Ellipsoids in
/// geodesy. Each ellipsoid comes with a list of fundamental geometric
/// characteristics (semi-major axis, flattening, name).
///
/// There are two ways a user can make use of the ellipsoids:
/// * 1. If the ellipsoid of choice is known at compile time, then the
///      template function/traits can be used. A (per-ellipsoid) specialized
///      traits class (i.e. ngpt::ellipsoid_traits) is used to implement the
///      basic properties of each ellipsoid.
///      E.g., it the ellipsoid of choice is GRS80, then:
///      double semi_major = ellipsoid_traits<ellipsoid::grs80>::a;
///      or
///      double N = N<ellipsoid::grs80>(lat);
///      gives the normal radius of curvature at latitude lat.
///
/// * 2. If the ellipsoid of choice is only kown at runtime, then all relevant
///      computations/constants can be accessed via the Ellipsoid class; i.e.
///      Ellipsoid e (ellipsoid::grs80);
///      double semi_major = e.semi_major();
///      double N = e.N(lat);
///
/// @example test_ellipsoid.cpp
///
/// [1] H. Moritz, GEODETIC REFERENCE SYSTEM 1980,
/// https://geodesy.geology.ohio-state.edu/course/refpapers/00740128.pdf

#ifndef __REFERENCE_ELLIPSOID__
#define __REFERENCE_ELLIPSOID__

#include <cmath>

namespace ngpt {

/// @brief core namespace holds the core of ellipsoid-related functions
namespace core {

// @brief Table of coefficients for meridian arc computation
double meridian_arc_length_impl1(double a, double f, double lat) noexcept;
double meridian_arc_length_impl2(double a, double f, double lat) noexcept;
double meridian_arc_length_impl3(double a, double f, double lat) noexcept;
double meridian_arc_length_impl4(double a, double f, double lat) noexcept;

/// @brief Compute the squared eccentricity.
///
/// Compute the squared (first) eccentricity (i.e. e^2) given the
/// flattening of an ellipsoid, aka \f$ e^2 = (2-f)*f \f$
/// @param[in] f flattening
/// @return      squared eccentricity
constexpr double eccentricity_squared(double f) noexcept {
  return (2e0 - f) * f;
}

/// @brief Compute the semi-minor axis (aka b)
///
/// Compute the semi-minor axis of an ellipsoid(i.e. e^2), given the
/// flattening and the semi-major axis, aka \f$ \beta = \alpha * (1-f) \f$
/// @param[in] f flattening
/// @param[in] a semi-major axis (meters)
/// @return      semi-minor axis (meters)
constexpr double semi_minor(double a, double f) noexcept {
  return a * (1e0 - f);
}

/// @brief the linear eccentricity
///
/// Compute the linear eccentricity of an ellipsoid (i.e. E) from the formula
/// \f$ E = \sqrt{a^2 - b^2} \f$
/// @param[in] f flattening
/// @param[in] a semi-major axis (meters)
/// @return      semi-minor axis (meters)
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
/// Compute the polar radius of curvature of an ellipsoid (i.e. E) from the
/// formula \f$ c = a^2 / b\f$
/// @param[in] f flattening
/// @param[in] a semi-major axis (meters)
/// @return      semi-minor axis (meters)
constexpr double polar_radius_of_curvature(double a, double f) noexcept {
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
} // namespace core

/// @brief A list of well-known reference ellipsoids.
///
/// For each of these reference ellipsoids, a series of traits (i.e. geometric
/// characteristics) will be specialized later on, using the template class
/// ngpt::ellipsoid_traits.
///
enum class ellipsoid : char { grs80, wgs84, pz90 };

/// @brief A class to hold ellipsoid traits (generic case).
///
/// A (template class) to hold specialized geometric quantities for each
/// of the eumerated elements (i.e. reference ellipsoids) in the
/// ngpt::ellipsoid enum.
/// I.e., to make any element of ngpt::ellipsoid usable, specialize this
/// (trait) class.
///
/// @tparam E  The reference ellipsoid to be specialized (i.e. one of
///            ngpt::ellipsoid).
template <ellipsoid E> struct ellipsoid_traits {};

/// @brief A class to hold traits for the GRS-80 (i.e. ngpt::ellispoid::grs80)
/// reference ellipsoid.
///
/// @see https://en.wikipedia.org/wiki/GRS_80
template <> struct ellipsoid_traits<ellipsoid::grs80> {
  /// Semi-major axis (m).
  static constexpr double a{6378137e0};
  /// Flattening.
  static constexpr double f{1e0 / 298.257222101e0};
  /// Reference ellipsoid name.
  static constexpr const char *n{"GRS80"};
};

/// @brief A class to hold traits for the WGS-84 (i.e. ngpt::ellispoid::wgs84)
/// reference ellipsoid.
///
/// @see https://en.wikipedia.org/wiki/World_Geodetic_System
template <> struct ellipsoid_traits<ellipsoid::wgs84> {
  /// Semi-major axis (m).
  static constexpr double a{6378137e0};
  /// Flattening.
  static constexpr double f{1e0 / 298.257223563e0};
  /// Reference ellipsoid name.
  static constexpr const char *n{"WGS84"};
};

/// @brief A class to hold traits for the PZ-90 (i.e. ngpt::ellispoid::pz90)
/// reference ellipsoid.
///
/// @see
/// http://www.navipedia.net/index.php/Reference_Frames_in_GNSS#GLONASS_reference_frame_PZ-90
template <> struct ellipsoid_traits<ellipsoid::pz90> {
  /// Semi-major axis (m).
  static constexpr double a{6378136e0};
  /// Flattening.
  static constexpr double f{1e0 / 298.257839303e0};
  /// Reference ellipsoid name.
  static constexpr const char *n{"PZ90"};
};

template <ellipsoid E>
double meridian_arc_length(double lat, int alg = 0) noexcept {
  constexpr double a = ellipsoid_traits<E>::a;
  constexpr double f = ellipsoid_traits<E>::f;
  switch (alg) {
  case 0:
    return core::meridian_arc_length_impl1(a, f, lat);
  case 1:
    return core::meridian_arc_length_impl2(a, f, lat);
  case 2:
    return core::meridian_arc_length_impl3(a, f, lat);
  default:
    return core::meridian_arc_length_impl1(a, f, lat);
  }
}

/// @brief Compute the squared eccentricity.
///
/// @tparam E  The reference ellipsoid (i.e. one of ngpt::ellipsoid).
/// @return    Eccentricity squared.
/// @throw     Does not throw.
///
/// @see ngpt::core::eccentricity_squared
template <ellipsoid E> constexpr double eccentricity_squared() noexcept {
  return core::eccentricity_squared(ellipsoid_traits<E>::f);
}

/// @brief Compute the linear eccentricity.
///
/// @tparam E  The reference ellipsoid (i.e. one of ngpt::ellipsoid).
/// @return    linear eccentricity
/// @throw     Does not throw.
///
/// @see ngpt::core::linear_eccentricity

template <ellipsoid E>
#if defined(__GNUC__) && !defined(__llvm__)
constexpr
#endif
    double
    linear_eccentricity() noexcept {
  return core::linear_eccentricity(ellipsoid_traits<E>::a,
                                   ellipsoid_traits<E>::f);
}

/// @brief Compute the semi-minor axis (b).
///
/// @tparam E The reference ellipsoid (i.e. one of ngpt::ellipsoid).
/// @return   Semi-minor axis of the reference ellipsoid in meters.
///
/// @see ngpt::core::semi_minor
template <ellipsoid E> constexpr double semi_minor() noexcept {
  return core::semi_minor(ellipsoid_traits<E>::a, ellipsoid_traits<E>::f);
}

/// @brief Compute the polar radius of curvature (c).
///
/// @tparam E The reference ellipsoid (i.e. one of ngpt::ellipsoid).
/// @return   Polar radius of curvature of the reference ellipsoid in meters.
///
/// @see ngpt::core::polar_radius_of_curvature
template <ellipsoid E> constexpr double polar_radius_of_curvature() noexcept {
  return core::polar_radius_of_curvature(ellipsoid_traits<E>::a,
                                         ellipsoid_traits<E>::f);
}

/// @brief Compute the normal radius of curvature at a given latitude (on a
///        reference ellipsoid).
///
/// @param[in] lat The latitude in radians.
/// @tparam    E   The reference ellipsoid (i.e. one of ngpt::ellipsoid).
/// @return        The normal radius of curvature in meters.
//
/// @see ngpt::core::N
template <ellipsoid E> double N(double lat) noexcept {
  return core::N(lat, ellipsoid_traits<E>::a, semi_minor<E>());
}

/// @brief Compute the meridional radii of curvature at a given latitude (on a
///        reference ellipsoid).
///
/// @param[in] lat The latitude in radians.
/// @tparam    E   The reference ellipsoid (i.e. one of ngpt::ellipsoid).
/// @return        The meridional radius of curvature in meters.
///
/// @see ngpt::core::M
template <ellipsoid E> double M(double lat) noexcept {
  return core::M(lat, ellipsoid_traits<E>::a, semi_minor<E>());
}

/// @brief Compute mean earth radius; in geophysics, the International Union of
///        Geodesy and Geophysics (IUGG) defines the mean radius (denoted R1)
///        to be: R1 = (2α + b) / 3
/// @tparam    E   The reference ellipsoid (i.e. one of ngpt::ellipsoid).
/// @return        The mean earth radius (as defined by IUGG for ellipsoid E)
///
/// @see https://en.wikipedia.org/wiki/Earth_radius#Mean_radius
template <ellipsoid E> constexpr double mean_earth_radius() noexcept {
  return 2e0 * ellipsoid_traits<E>::a / 3e0 + semi_minor<E>() / 3e0;
}

/// @brief Arc length of an infinitesimal element of the meridian
///
/// @param[in] lat  Latitude of the reference point in radians
/// @param[in] dlat Lattitude difference, i.e. arc length on the meridian
/// @return         Arc length (on meridian) in meters
/// @warning This formula is valid for infinitesimal latitude differences
/// @see https://en.wikipedia.org/wiki/Meridian_arc
template <ellipsoid E, typename T,
          typename = std::enable_if_t<std::is_floating_point<T>::value>>
T infinitesimal_meridian_arc(T lat, T dlat) noexcept {
  T M_ = core::M(lat, ellipsoid_traits<E>::a, semi_minor<E>());
  return M_ * dlat;
}

/// @brief Arc length on parallel
///
/// @param[in] lat  Latitude of the parallel (radians)
/// @param[in] dlon Longtitude difference (radians)
/// @return         Arc length (on parallel) in meters
template <ellipsoid E, typename T,
          typename = std::enable_if_t<std::is_floating_point<T>::value>>
T parallel_arc_length(T lat, T dlon) noexcept {
  const T N_ = core::N(lat, ellipsoid_traits<E>::a, semi_minor<E>());
  const T cosf = std::cos(lat);
  return N_ * cosf * dlon;
}

/// @class Ellipsoid
///
/// A class to represent a reference ellipsoid. An ellipsoid is defined by
/// two parameters, namely:
/// * semi-major axis, \f$ \alpha \f$ aka the equatorial radius of the ellipsoid
/// * flattening, f aka \f$ f = \frac{\alpha - \beta}{\alpha} \f$
///
/// Users can construct the commonly used ellipsoids in Geodesy (grs80,
/// wgs84 and pz90) via the ngpt::ellipsoid enums, or any other ellipsoid
/// of choice, by passing in the fundamental arguments (a and f).
class Ellipsoid {
public:
  /// @brief  Constructor from an ngpt::ellipsoid enum
  /// @param[in] e An ngpt::ellipsoid; fundamental geometric constants are
  ///              automatically assigned via the ngpt::ellipsoid_traits
  ///              class.
  explicit constexpr Ellipsoid(ellipsoid e) noexcept : __a(0e0), __f(0e0) {
    switch (e) {
    case ellipsoid::grs80:
      __a = ellipsoid_traits<ellipsoid::grs80>::a;
      __f = ellipsoid_traits<ellipsoid::grs80>::f;
      break;
    case ellipsoid::wgs84:
      __a = ellipsoid_traits<ellipsoid::wgs84>::a;
      __f = ellipsoid_traits<ellipsoid::wgs84>::f;
      break;
    case ellipsoid::pz90:
      __a = ellipsoid_traits<ellipsoid::pz90>::a;
      __f = ellipsoid_traits<ellipsoid::pz90>::f;
      break;
    }
  }

  /// @brief  User-defined instance.
  /// @param[in] a The semi-major axis (meters)
  /// @param[in] f The flattening
  constexpr Ellipsoid(double a, double f) noexcept : __a(a), __f(f){};

  /// @brief  Get the semi-major axis \f$ \alpha \f$
  /// @return The semi-major axis (meters)
  constexpr double semi_major() const noexcept { return __a; }

  /// @brief  Get the flattening
  /// @return The flattening
  constexpr double flattening() const noexcept { return __f; }

  /// @brief  Get the squared eccentricity \f$ e^2 \f$
  /// @return Eccentricity aquared
  /// @see ngpt::core::eccentricity_squared
  constexpr double eccentricity_squared() const noexcept {
    return core::eccentricity_squared(__f);
  }

  /// @brief  Get the semi-minor axis \f$ \beta \f$
  /// @return Semi-minor axis (meters)
  /// @see ngpt::core::semi_minor
  constexpr double semi_minor() const noexcept {
    return core::semi_minor(__a, __f);
  }

  /// @brief  Compute the normal radius of curvature at a given latitude
  ///
  /// @param[in] lat The latitude in radians.
  /// @return        The normal radius of curvature in meters.
  ///
  /// @see ngpt::core::N
  double N(double lat) const noexcept {
    return core::N(lat, __a, this->semi_minor());
  }

  /// @brief  Compute the meridional radii of curvature at a given latitude
  ///
  /// @param[in] lat The latitude in radians.
  /// @return        The meridional radius of curvature in meters.
  ///
  /// @see ngpt::core::M
  double M(double lat) const noexcept {
    return core::M(lat, __a, this->semi_minor());
  }

  /// @brief Compute mean earth radius; in geophysics, the International Union
  /// of
  ///        Geodesy and Geophysics (IUGG) defines the mean radius (denoted R1)
  ///        to be: R1 = (2α + b) / 3
  /// @tparam    E   The reference ellipsoid (i.e. one of ngpt::ellipsoid).
  /// @return        The mean earth radius (as defined by IUGG for ellipsoid E)
  ///
  /// @see https://en.wikipedia.org/wiki/Earth_radius#Mean_radius
  constexpr double mean_earth_radius() const noexcept {
    return 2e0 * __a / 3e0 + this->semi_minor() / 3e0;
  }

private:
  double __a, __f;
}; // class Ellipsoid

} // end namespace ngpt

#endif
