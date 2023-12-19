/** @file units.hpp
 * A list of frequently used geodetic functions for unit conversion, mostly 
 * targeting angular units.
 */

#ifndef __DSO_GEODESY_UNITS_HPP__
#define __DSO_GEODESY_UNITS_HPP__

#include "geoconst.hpp"
#include <cmath>

namespace dso {

/** @brief Convert degrees to radians.
 * @param[in] degrees Angle in decimal [deg]
 * @return The (input) angle [rad]
 */
inline constexpr double deg2rad(double degrees) noexcept {
  constexpr const double f = detail::AngleUnitTraits<
      detail::AngleUnit::Degrees>::to_units<detail::AngleUnit::Radians>();
  return degrees * f;
}

/** @brief Convert radians to degrees.
 *  @param[in] radians Angle in [rad]
 *  @return The (input) angle in decimal [deg]
 */
inline constexpr double rad2deg(double radians) noexcept {
  constexpr const double f = detail::AngleUnitTraits<
      detail::AngleUnit::Radians>::to_units<detail::AngleUnit::Degrees>();
  return radians * f;
}

/** @brief Convert radians to seconds (of degree).
 *  @param[in] radians Angle in [rad]
 *  @return The (input) angle in [sec] (of degrees)
 */
inline constexpr double rad2sec(double radians) noexcept {
  constexpr const double f = detail::AngleUnitTraits<
      detail::AngleUnit::Radians>::to_units<detail::AngleUnit::Seconds>();
  return radians * f;
}

/** @brief Convert seconds (of degree) to radians.
 *  @param[in] seconds Angle in [sec] of degree
 *  @return The (input) angle in [rad]
 */
inline constexpr double sec2rad(double seconds) noexcept {
  constexpr const double f = detail::AngleUnitTraits<
      detail::AngleUnit::Seconds>::to_units<detail::AngleUnit::Radians>();
  return seconds * f;
}

/** Normalize angle in the range [0, 2π]/[0,360]
 * 
 * @tparam U Units of input angle, AngleUnits::Radians, or AngleUnits::Degrees
 *         or ...
 * @param[in] a Angle in units of U
 * @return Normalized angle in the range 0 to 1-cycle, in units of U
 */
template <detail::AngleUnit U = detail::AngleUnit::Radians>
inline double norm_angle(double a) noexcept {
  constexpr const double circle = detail::AngleUnitTraits<U>::full_circle();
  a = std::fmod(a, circle);
  const double r[] = {a, a + circle};
  return r[a < 0e0];
}

template <detail::AngleUnit U = detail::AngleUnit::Radians>
inline double anp(double a) noexcept { return norm_angle<U>(a); }

/** Normalize angle in the range [-1/2 to 1/2) of circle.
 *
 * @tparam U Units of input angle, AngleUnits::Radians, or AngleUnits::Degrees
 *         or ...
 * @param[in] a Angle in units of U
 * @return Normalized angle in the range -1/2 to 1/2 -cycle, in units of U
 */
template <detail::AngleUnit U = detail::AngleUnit::Radians>
inline double anpm(double a) noexcept {
  constexpr const double circle = detail::AngleUnitTraits<U>::full_circle();
  constexpr const double halfCircle =
      detail::AngleUnitTraits<U>::full_circle() / 2e0;
  a = std::fmod(a, circle);
  const double r[] = {a, a - std::copysign(circle, a)};
  return r[std::abs(a) >= halfCircle];
}

/** @brief Decimal to hexicondal degrees.
 *
 * @param[in]  decimal_deg The decimal degrees.
 * @param[out] deg         Integer degrees.
 * @param[out] min         Integer minutes.
 * @param[out] sec         The fractional seconds.
 * @param[out] sign        Either -1 or +1; this represents the sign of the
 *                         hexicondal degrees
 *
 * @note
 * Why use a seperate variable for sign? Well, if we didn't, we could e.g. 
 * make deg positive or negative depending on decimal_deg. BUT, if deg is 
 * zero, then we have a problem cause -0 cannot be (usually) represented, so 
 * the signs of decimal_deg and deg would not be the same! Hence, deg, min 
 * and sec are always >= 0 and the sign of the angle is the sign of sign 
 * (variable).
 *
 *  E.g.
 *  dso::decd2hexd(10e0, deg1, min1, sec1, sgn1);
 *  dso::decd2hexd(-10e0, deg2, min2, sec2, sgn2);
 *  assert(deg1==deg2 && (min1==min2 && sec1==sec2));
 *  assert(sgn1==-sgn2);
 */
inline constexpr void decd2hexd(double decimal_deg, int &deg, int &min,
                                double &sec, int &sign) noexcept {
  const double decdeg = std::abs(decimal_deg);
  deg = static_cast<int>(decdeg);
  min = static_cast<int>((decdeg - static_cast<double>(deg)) *
                         static_cast<double>(60e0));
  sec = decdeg - (static_cast<double>(deg) + static_cast<double>(min) / 60e0);
  sec *= 3600e0;

  sign = static_cast<int>(std::copysign(1e0, decimal_deg));
  return;
}

/** @brief Hexicondal degrees to decimal degrees.
 *
 * @param[in] deg  Integer degrees.
 * @param[in] min  Integer minutes.
 * @param[in] sec  The fractional seconds.
 * @param[in] sign The sign of the hexicondal degrees; that is any integer with
 *                 the correct sign (we only consider the sign of the parameter
 *                 not its value).
 * @return         The angle in decimal degrees.
 * @throw          Does not throw
 *
 * @note If the angle is negative, only the sign parameter should be negative;
 *  (deg parameter could also be negative; the function will only use its
 *  absolute value, disregarding the sign).
 *  if so, the (decimal) degrees returned will also be negative. The
 *  parameters min and sec are not checked for their sign, they should
 *  *ALWAYS* be positive.
 *  Why do we need the sign parameter? Well, if deg is zero, but the
 *  degrees are negative (e.g. -0deg 10min 10.10sec), how could we know
 *  the sign? There is no -0; so we need to mark the sign via the sign
 *  parameter.
 *
 * E.g.
 *  dso::decd2hexd(10e0, deg1, min1, sec1, sgn1);
 *  dso::decd2hexd(-10e0, deg2, min2, sec2, sgn2);
 *  assert(deg1==deg2 && (min1==min2 && sec1==sec2));
 *  assert(sgn1==-sgn2);
 *  a1 = dso::hexd2decd(deg1, min1, sec1, sgn1);
 *  a2 = dso::hexd2decd(deg2, min2, sec2, sgn2);
 *  assert(a1 == -a2);
 *  // note the sign of deg parameter is not considered!
 *  a2 = dso::hexd2decd(-deg1, min1, sec1, sgn1);
 *  assert(a1==a2);
 */
inline constexpr double hexd2decd(int deg, int min, double sec,
                                  int sign = 1) noexcept {
  double angle = static_cast<double>(std::abs(deg)) +
                 (static_cast<double>(min) + sec / 60e0) / 60e0;
  return std::copysign(angle, (double)sign);
}

/** @brief Hexicondal degrees to radians.
 *
 * @param[in] deg  Integer degrees.
 * @param[in] min  Integer minutes.
 * @param[in] sec  The fractional seconds.
 * @param[in] sign The sign of the hexicondal degrees; that is any integer with
 *                 the correct sign (we only consider the sign of the parameter
 *                 not its value).
 * @return         The angle in [rad].
 *
 * @note If the angle is negative, only the sign parameter should be negative;
 *  (deg parameter could also be negative; the function will only use its
 *  absolute value, disregarding the sign).
 *  if so the (decimal) degrees returned will also be negative. The parameters 
 *  min and sec are not checked for their sign, they should *ALWAYS* be 
 *  positive.
 *  Why do we need the sign parameter? Well, if deg is zero, but the degrees 
 *  are negative (e.g. -0deg 10min 10.10sec), how could we know the sign? 
 *  There is no -0; so we need to mark the sign via the sign parameter.
 *
 * E.g.
 *  a1 = dso::hexd2rad(deg1, min1, sec1, sgn1);
 *  a2 = dso::hexd2rad(deg2, min2, sec2, sgn2);
 *  assert(a1==-a2);
 *  // note the sign of deg parameter is not considered!
 *  a1 = dso::hexd2rad(-deg1, min1, sec1, sgn1);
 *  assert(a1==-a2);
 */
inline constexpr double hexd2rad(int deg, int min, double sec,
                                 int sign = 1) noexcept {
  return deg2rad(hexd2decd(deg, min, sec, sign));
}

/** @brief Radians to hexicondal degrees.
 *
 * @param[in]  radians     An angle in [rad].
 * @param[out] deg         Integer degrees.
 * @param[out] min         Integer minutes.
 * @param[out] sec         The fractional seconds.
 * @param[out] sign        Either -1 or +1; this represents the sign of the
 *                         hexicondal degrees
 *
 * @note In case a negative angle is given, then the (output) degrees are also
 *  going to be negative.
 *  Why use a seperate variable for sign? Well, if we didn't, we could e.g. 
 *  make deg positive or negative depending on decimal_deg. BUT, if deg is 
 *  zero, then we have a problem cause -0 cannot be (usually) represented, so 
 *  the signs of decimal_deg and deg would not be the same!
 */
inline constexpr void rad2hexd(double radians, int &deg, int &min, double &sec,
                               int &sign) noexcept {
  return decd2hexd(rad2deg(radians), deg, min, sec, sign);
}

} /* namespace dso */

#endif
