///
/// @file units.hpp
///
/// @brief A list of frequently used geodetic functions for unit conversion, 
///        mostly targeted on angular units.
///

#ifndef __NGPT_GEODESY_UNITS_HPP__
#define __NGPT_GEODESY_UNITS_HPP__

#include "geoconst.hpp"
#include <cassert>
#include <cmath>

namespace ngpt {

/// @brief Convert degrees to radians.
/// @tparam    T       Any floating type
/// @param[in] degrees Angle in decimal degrees
/// @return            The (input) angle in radians
template <typename T> constexpr T deg2rad(T degrees) noexcept {
  return degrees * DEG2RAD;
}

/// @brief Convert radians to degrees.
/// @tparam    T       Any floating type
/// @param[in] radians Angle in radians
/// @return            The (input) angle in decimal degrees
template <typename T> constexpr T rad2deg(T radians) noexcept {
  return radians * RAD2DEG;
}

/// @brief Convert radians to seconds (of degree).
/// @tparam    T       Any floating type
/// @param[in] radians Angle in radians
/// @return            The (input) angle in seconds (of degrees)
template <typename T> constexpr T rad2sec(T radians) noexcept {
  return (radians * RAD2DEG) * 3600e0;
}

/// @brief Normalize angle.
///
/// Normalize an angle in the interval [lower, upper).
///
/// @tparam    T      Any floating point type for input and results.
/// @param[in] angle  The angle to normalize (note that the unit should be the
///                   same as in lower and upper parameters).
/// @param[in] lower  lower bound (inclusive). Default value is 0
/// @param[in] upper  upper bound (exclusive). Default is 2* \f$\pi\f$
///
/// @note  It is not always needed to use this function to normalize an angle.
///        If i.e. the angle is a result of a function that returns values in
///        the range (-\f$\pi\f$ , +\f$\pi\f$] radians, and we need to normalize
///        in the range [0, 2\f$\pi\f$), then we can use: angle =
///        std::fmod(angle+2\f$\pi\f$, 2\f$\pi\f$).
template <typename T,
          typename = std::enable_if_t<std::is_floating_point<T>::value>>
T normalize_angle(T angle, T lower = 0e0, T upper = D2PI) noexcept {
  assert(lower < upper);

  // std::fmod can do the job for a range with lower limit >= 0 (e.g. 0 to 2\pi)
  // but will fail for e.g. -\pi to \pi
  if (lower >= 0e0)
    return std::fmod(angle + upper, upper);

  double res{angle};
  if (angle > upper || angle == lower)
    angle = lower + std::fmod(std::abs(angle + upper),
                              std::abs(lower) + std::abs(upper));
  if (angle < lower || angle == upper)
    angle = upper - std::fmod(std::abs(angle - lower),
                              std::abs(lower) + std::abs(upper));

  res = (res == upper) ? (lower) : (angle);
  return res;
}

/// @brief Decimal to hexicondal degrees.
///
/// @tparam     T           Any floating point type for input and results.
/// @param[in]  decimal_deg The decimal degrees.
/// @param[out] deg         Integer degrees.
/// @param[out] min         Integer minutes.
/// @param[out] sec         The fractional seconds.
/// @param[out] sign        Either -1 or +1; this represents the sign of the
///                         hexicondal degrees
/// @throw                  Does not throw
///
/// @note
///       Why use a seperate variable for sign? Well, if we didn't, we could
///       e.g. make deg positive or negative depending on decimal_deg. BUT, if
///       deg is zero, then we have a problem cause -0 cannot be (usually)
///       represented, so the signs of decimal_deg and deg would not be the
///       same! Hence, deg, min and sec are always >= 0 and the sign of the
///       angle is the sign of sign (variable).
///  E.g.
///  ngpt::decd2hexd(10e0, deg1, min1, sec1, sgn1);
///  ngpt::decd2hexd(-10e0, deg2, min2, sec2, sgn2);
///  assert(deg1==deg2 && (min1==min2 && sec1==sec2));
///  assert(sgn1==-sgn2);
template <typename T,
          typename = std::enable_if_t<std::is_floating_point<T>::value>>
constexpr void decd2hexd(T decimal_deg, int &deg, int &min, T &sec,
                         int &sign) noexcept {
  T decdeg{std::abs(decimal_deg)};
  deg = static_cast<int>(decdeg);
  min = static_cast<int>((decdeg - static_cast<T>(deg)) * static_cast<T>(60e0));
  sec = decdeg - (static_cast<T>(deg) + static_cast<T>(min) / 60e0);
  sec *= 3600e0;

  sign = static_cast<int>(std::copysign(1e0, decimal_deg));
  return;
}

/// @brief Hexicondal degrees to decimal degrees.
///
/// @tparam     T   Any floating point type for input and results.
/// @param[in] deg  Integer degrees.
/// @param[in] min  Integer minutes.
/// @param[in] sec  The fractional seconds.
/// @param[in] sign The sign of the hexicondal degrees; that is any integer with
///                 the correct sign (we only consider the sign of the parameter
///                 not its value).
/// @return         The angle in decimal degrees.
/// @throw          Does not throw
///
/// @note If the angle is negative, only the sign parameter should be negative;
///       (deg parameter could also be negative; the function will only use its
///       absolute value, disregarding the sign).
///       if so, the (decimal) degrees returned will also be negative. The
///       parameters min and sec are not checked for their sign, they should
///       *ALWAYS* be positive.
///       Why do we need the sign parameter? Well, if deg is zero, but the
///       degrees are negative (e.g. -0deg 10min 10.10sec), how could we know
///       the sign? There is no -0; so we need to mark the sign via the sign
///       parameter.
/// E.g.
///  ngpt::decd2hexd(10e0, deg1, min1, sec1, sgn1);
///  ngpt::decd2hexd(-10e0, deg2, min2, sec2, sgn2);
///  assert(deg1==deg2 && (min1==min2 && sec1==sec2));
///  assert(sgn1==-sgn2);
///  a1 = ngpt::hexd2decd(deg1, min1, sec1, sgn1);
///  a2 = ngpt::hexd2decd(deg2, min2, sec2, sgn2);
///  assert(a1 == -a2);
///  // note the sign of deg parameter is not considered!
///  a2 = ngpt::hexd2decd(-deg1, min1, sec1, sgn1);
///  assert(a1==a2);
template <typename T,
          typename = std::enable_if_t<std::is_floating_point<T>::value>>
constexpr T hexd2decd(int deg, int min, T sec, int sign = 1) noexcept {
  T angle{static_cast<T>(std::abs(deg)) +
          (static_cast<T>(min) + sec / 60e0) / 60e0};
  return std::copysign(angle, (T)sign);
}

/// @brief Hexicondal degrees to radians.
///
/// @tparam     T   Any floating point type for input and results.
/// @param[in] deg  Integer degrees.
/// @param[in] min  Integer minutes.
/// @param[in] sec  The fractional seconds.
/// @param[in] sign The sign of the hexicondal degrees; that is any integer with
///                 the correct sign (we only consider the sign of the parameter
///                 not its value).
/// @return         The angle in radians.
/// @throw          Does not throw
///
/// @note If the angle is negative, only the sign parameter should be negative;
///       (deg parameter could also be negative; the function will only use its
///       absolute value, disregarding the sign).
///       if so the (decimal) degrees returned will also be negative. The
///       parameters min and sec are not checked for their sign, they should
///       *ALWAYS* be positive.
///       Why do we need the sign parameter? Well, if deg is zero, but the
///       degrees are negative (e.g. -0deg 10min 10.10sec), how could we know
///       the sign? There is no -0; so we need to mark the sign via the sign
///       parameter.
/// E.g.
///  a1 = ngpt::hexd2rad(deg1, min1, sec1, sgn1);
///  a2 = ngpt::hexd2rad(deg2, min2, sec2, sgn2);
///  assert(a1==-a2);
///  // note the sign of deg parameter is not considered!
///  a1 = ngpt::hexd2rad(-deg1, min1, sec1, sgn1);
///  assert(a1==-a2);
template <typename T,
          typename = std::enable_if_t<std::is_floating_point<T>::value>>
constexpr T hexd2rad(int deg, int min, T sec, int sign = 1) noexcept {
  return deg2rad(hexd2decd(deg, min, sec, sign));
}

/// @brief Radians to hexicondal degrees.
///
/// @tparam     T           Any floating point type for input and results.
/// @param[in]  radians     An angle in radians.
/// @param[out] deg         Integer degrees.
/// @param[out] min         Integer minutes.
/// @param[out] sec         The fractional seconds.
/// @param[out] sign        Either -1 or +1; this represents the sign of the
///                         hexicondal degrees
/// @throw                  Does not throw
///
/// @note In case a negative angle is given, then the (output) degrees are also
///       going to be negative.
///       Why use a seperate variable for sign? Well, if we didn't, we could
///       e.g. make deg positive or negative depending on decimal_deg. BUT, if
///       deg is zero, then we have a problem cause -0 cannot be (usually)
///       represented, so the signs of decimal_deg and deg would not be the
///       same!
template <typename T,
          typename = std::enable_if_t<std::is_floating_point<T>::value>>
constexpr void rad2hexd(T radians, int &deg, int &min, T &sec,
                        int &sign) noexcept {
  return decd2hexd(rad2deg(radians), deg, min, sec, sign);
}

}// ngpt

#endif
