///
/// @file geoconst.hpp
///
/// @brief Fundamental constants frequently used within geodetic calculations.
///

#ifndef __NGPT_GEOCONST_HPP__
#define __NGPT_GEOCONST_HPP__
#include <cmath>

namespace ngpt
{
  /// The value of pi.
  constexpr double DPI  { std::atan(1e0)*4e0 };
  
  /// The value of 2 * pi.
  constexpr double D2PI { 2e0 * DPI };

  /// Degrees to Radians coefficient.
  constexpr double DEG2RAD { DPI / 180e0 };
  
  /// Radians to Degrees coefficient.
  constexpr double RAD2DEG { 180e0 / DPI };

} // namespace ngpt

#endif
