/** @file
 * Fundamental constants frequently used within geodetic calculations.
 */

#ifndef __DSO_GEODESY_CONSTS_HPP__
#define __DSO_GEODESY_CONSTS_HPP__

#include <cmath>

namespace dso {

namespace detail {
/** Type of angular units */
enum class AngleUnit : char { Radians, Degrees, Seconds };

/** Traits of default angular units (AngleUnit) i.e. AngleUnit::Radians */
template <AngleUnit U> struct AngleUnitTraits {
  /** A full circle in given angular units (default=Radians) */
  static constexpr double full_circle() { return 2e0 * M_PI; }
  /** half a circle in given angular units (default=Radians) */
  static constexpr double half_circle() { return M_PI; }
  /** Factor to convert to some other angular unit */
  template<AngleUnit T> static constexpr double to_units() noexcept {
    return AngleUnitTraits<T>::half_circle() / half_circle();
  }
};

/** Traits of degrees angular units (AngleUnit) i.e. AngleUnit::Degrees */
template <> struct AngleUnitTraits<AngleUnit::Degrees> {
  /** A full circle in given angular units, i.e. Degrees */
  static constexpr double full_circle() { return 360e0; }
  /** Half a circle in given angular units, i.e. Degrees */
  static constexpr double half_circle() { return 180e0; }
  /** Factor to convert to some other angular unit */
  template<AngleUnit T> static constexpr double to_units() noexcept {
    return AngleUnitTraits<T>::half_circle() / half_circle();
  }
};

/** Traits of degrees angular units (AngleUnit) i.e. AngleUnit::Seconds */
template <> struct AngleUnitTraits<AngleUnit::Seconds> {
  /** A full circle in given angular units, i.e. Seconds (of degree)*/
  static constexpr double full_circle() { return (double)(360L * 60L * 60L); }
  /** Half a circle in given angular units, i.e. Seconds (of degree) */
  static constexpr double half_circle() { return (double)(180L * 60L * 60L); }
  /** Factor to convert to some other angular unit */
  template<AngleUnit T> static constexpr double to_units() noexcept {
    return AngleUnitTraits<T>::half_circle() / half_circle();
  }
};
} /* namespace detail */

/** The value of \f$\pi\f$. */
#if defined(__GNUC__) && !defined(__llvm__)
constexpr const double DPI = std::atan(1e0) * 4e0;
#else
constexpr const double DPI = M_PI;
#endif

/** The value of \f$2*\pi\f$ */
constexpr const double D2PI =
    detail::AngleUnitTraits<detail::AngleUnit::Radians>::full_circle();

/** Degrees to Radians coefficient. */
constexpr const double DEG2RAD = detail::AngleUnitTraits<
    detail::AngleUnit::Degrees>::to_units<detail::AngleUnit::Radians>();

/** Radians to Degrees coefficient. */
constexpr const double RAD2DEG = detail::AngleUnitTraits<
    detail::AngleUnit::Radians>::to_units<detail::AngleUnit::Degrees>();

// TODO ??
constexpr const double TURNAS = 1296000e0;

} /* namespace dso */

#endif
