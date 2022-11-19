/// @file ellipsoid.hpp
/// @brief Definition and basic setup of an ellipsoid (kinda) class.
///
/// This file defines a list of frequently used (referennce) Ellipsoids in
/// geodesy. Each ellipsoid comes with a list of fundamental geometric
/// characteristics (semi-major axis, flattening, name).
///
/// There are two ways a user can make use of the ellipsoids:
/// * 1. If the ellipsoid of choice is known at compile time, then the
///      template function/traits can be used. A (per-ellipsoid) specialized
///      traits class (i.e. dso::ellipsoid_traits) is used to implement the
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
/// Note that the semi-major axis is sometimes reffered to as
/// "equatorial radius" and the semi-minor axis is sometimes reffered to as
/// "polar radius".
/// In the following and when no other qualification is used, latitude is used
/// for Geodetic latitude.
///
/// @example test_ellipsoid.cpp
///
/// [1] H. Moritz, GEODETIC REFERENCE SYSTEM 1980,
/// https://geodesy.geology.ohio-state.edu/course/refpapers/00740128.pdf
///
/// [2] Charles F. F. Karney, Algorithms for geodesics, J Geod (2013) 87:43–55
/// [3] https://en.wikipedia.org/wiki/Latitude

#ifndef __REFERENCE_ELLIPSOID_HPP__
#define __REFERENCE_ELLIPSOID_HPP__

#include "ellipsoid_core.hpp"
#include <type_traits>

namespace dso {

/// @brief A list of well-known reference ellipsoids.
///
/// For each of these reference ellipsoids, a series of traits (i.e. geometric
/// characteristics) will be specialized later on, using the template class
/// dso::ellipsoid_traits.
///
enum class ellipsoid : char { grs80, wgs84, pz90 };

/// @brief A class to hold ellipsoid traits (generic case).
///
/// A (template class) to hold specialized geometric quantities for each
/// of the eumerated elements (i.e. reference ellipsoids) in the
/// dso::ellipsoid enum.
/// I.e., to make any element of dso::ellipsoid usable, specialize this
/// (trait) class.
///
/// @tparam E  The reference ellipsoid to be specialized (i.e. one of
///            dso::ellipsoid).
template <ellipsoid E> struct ellipsoid_traits {};

/// @brief A class to hold traits for the GRS-80 (i.e. dso::ellispoid::grs80)
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

/// @brief A class to hold traits for the WGS-84 (i.e. dso::ellispoid::wgs84)
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

/// @brief A class to hold traits for the PZ-90 (i.e. dso::ellispoid::pz90)
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

/// @brief Compute the squared eccentricity.
///
/// @tparam E  The reference ellipsoid (i.e. one of dso::ellipsoid).
/// @return    Eccentricity squared.
/// @throw     Does not throw.
///
/// @see dso::core::eccentricity_squared
template <ellipsoid E> constexpr double eccentricity_squared() noexcept {
  return core::eccentricity_squared(ellipsoid_traits<E>::f);
}

/// @brief Compute the linear eccentricity.
///
/// @tparam E  The reference ellipsoid (i.e. one of dso::ellipsoid).
/// @return    linear eccentricity
/// @throw     Does not throw.
///
/// @see dso::core::linear_eccentricity
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
/// @tparam E The reference ellipsoid (i.e. one of dso::ellipsoid).
/// @return   Semi-minor axis of the reference ellipsoid in meters.
///
/// @see dso::core::semi_minor
template <ellipsoid E> constexpr double semi_minor() noexcept {
  return core::semi_minor(ellipsoid_traits<E>::a, ellipsoid_traits<E>::f);
}

/// @brief Compute the polar radius of curvature (c).
///
/// @tparam E The reference ellipsoid (i.e. one of dso::ellipsoid).
/// @return   Polar radius of curvature of the reference ellipsoid in meters.
///
/// @see dso::core::polar_radius_of_curvature
template <ellipsoid E> constexpr double polar_radius_of_curvature() noexcept {
  return core::polar_radius_of_curvature(ellipsoid_traits<E>::a,
                                         ellipsoid_traits<E>::f);
}

/// @brief Compute the third flattening
///
/// @tparam E The reference ellipsoid (i.e. one of dso::ellipsoid).
/// @return   The third flattening
///
/// @see dso::core::third_flattening
template <ellipsoid E> constexpr double third_flattening() noexcept {
  return core::third_flattening(ellipsoid_traits<E>::f);
}

/// @brief Compute the geocentric latitude at some geodetic latitude on the
///        ellipsoid
///
/// @tparam    E   The reference ellipsoid (i.e. one of dso::ellipsoid).
/// @param[in] lat The geodetic latitude [rad]
/// @return        The geocentric latitude [rad]
/// @throw     Does not throw.
///
/// @see dso::core::geocentric_latitude
template <ellipsoid E>
#if defined(__GNUC__) && !defined(__llvm__)
constexpr
#endif
    double
    geocentric_latitude(double lat) noexcept {
  return core::geocentric_latitude(lat, ellipsoid_traits<E>::f);
}

/// @brief Compute the parametric or reduced latitude at some geodetic latitude
///        on the ellipsoid
///
/// @tparam    E   The reference ellipsoid (i.e. one of dso::ellipsoid).
/// @param[in] lat The geodetic latitude in radians
/// @return        The parametric or reduced latitude in radians
/// @throw     Does not throw.
///
/// @see dso::core::reduced_latitude
template <ellipsoid E>
#if defined(__GNUC__) && !defined(__llvm__)
constexpr
#endif
    double
    reduced_latitude(double lat) noexcept {
  return core::reduced_latitude(lat, ellipsoid_traits<E>::f);
}

/// @brief Compute the normal radius of curvature at a given latitude (on a
///        reference ellipsoid).
///
/// @param[in] lat The latitude in radians.
/// @tparam    E   The reference ellipsoid (i.e. one of dso::ellipsoid).
/// @return        The normal radius of curvature in meters.
//
/// @see dso::core::N
template <ellipsoid E> double N(double lat) noexcept {
  return core::N(lat, ellipsoid_traits<E>::a, semi_minor<E>());
}

/// @brief Compute the meridional radii of curvature at a given latitude (on a
///        reference ellipsoid).
///
/// @param[in] lat The latitude in radians.
/// @tparam    E   The reference ellipsoid (i.e. one of dso::ellipsoid).
/// @return        The meridional radius of curvature in meters.
///
/// @see dso::core::M
template <ellipsoid E> double M(double lat) noexcept {
  return core::M(lat, ellipsoid_traits<E>::a, semi_minor<E>());
}

/// @brief Compute mean earth radius; in geophysics, the International Union of
///        Geodesy and Geophysics (IUGG) defines the mean radius (denoted R1)
///        to be: \f$ R1 = (2\alpha + b) / 3 \f$
/// @tparam    E   The reference ellipsoid (i.e. one of dso::ellipsoid).
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
/// wgs84 and pz90) via the dso::ellipsoid enums, or any other ellipsoid
/// of choice, by passing in the fundamental arguments (a and f).
class Ellipsoid {
public:
  /// @brief  Constructor from an dso::ellipsoid enum
  /// @param[in] e An dso::ellipsoid; fundamental geometric constants are
  ///              automatically assigned via the dso::ellipsoid_traits
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
  /// @see dso::core::eccentricity_squared
  constexpr double eccentricity_squared() const noexcept {
    return core::eccentricity_squared(__f);
  }

  /// @brief  Get the semi-minor axis \f$ \beta \f$
  /// @return Semi-minor axis (meters)
  /// @see dso::core::semi_minor
  constexpr double semi_minor() const noexcept {
    return core::semi_minor(__a, __f);
  }

  /// @brief Get the third flattening \f$ n \f$
  /// @return The third flattening
  /// @see dso::core::third_flattening
  constexpr double third_flattening() const noexcept {
    return core::third_flattening(__f);
  }

  /// @brief Compute the geocentric latitude at some (geodetic) latitude
  /// @param[in] lat The geodecit latitude
  /// @return        The geocentric latitude at lat
  double geocentric_latitude(double lat) const noexcept {
    return core::geocentric_latitude(lat, __f);
  }

  /// @brief Compute the reduced latitude at some (geodetic) latitude
  /// @param[in] lat The geodecit latitude
  /// @return        The reduced latitude at lat
  double reduced_latitude(double lat) const noexcept {
    return core::reduced_latitude(lat, __f);
  }

  /// @brief  Compute the normal radius of curvature at a given latitude
  ///
  /// @param[in] lat The latitude in radians.
  /// @return        The normal radius of curvature in meters.
  ///
  /// @see dso::core::N
  double N(double lat) const noexcept {
    return core::N(lat, __a, this->semi_minor());
  }

  /// @brief  Compute the meridional radii of curvature at a given latitude
  ///
  /// @param[in] lat The latitude in radians.
  /// @return        The meridional radius of curvature in meters.
  ///
  /// @see dso::core::M
  double M(double lat) const noexcept {
    return core::M(lat, __a, this->semi_minor());
  }

  /// @brief Compute mean earth radius; in geophysics, the International Union
  /// of
  ///        Geodesy and Geophysics (IUGG) defines the mean radius (denoted R1)
  ///        to be: \f$ R1 = (2\alpha + b) / 3 \f$
  /// @tparam    E   The reference ellipsoid (i.e. one of dso::ellipsoid).
  /// @return        The mean earth radius (as defined by IUGG for ellipsoid E)
  ///
  /// @see https://en.wikipedia.org/wiki/Earth_radius#Mean_radius
  constexpr double mean_earth_radius() const noexcept {
    return 2e0 * __a / 3e0 + this->semi_minor() / 3e0;
  }

private:
  double __a, __f;
}; // Ellipsoid

} // namespace dso

#endif
