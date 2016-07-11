///
/// \file geoconst.hpp
///
/// \brief Fundamental constants frequently used within geodetic calculations.
///

#ifndef __NGPT_GEOCONST_HPP__
#define __NGPT_GEOCONST_HPP__

namespace ngpt
{
    /// The value of pi.
    constexpr double DPI  { 3.141592653589793238463 };
    
    /// The value of 2 * pi.
    constexpr double D2PI { 2 * DPI };

    /// Degrees to Radians coefficient.
    constexpr double DEG2RAD { DPI / 180.0e0 };
    
    /// Radians to Degrees coefficient.
    constexpr double RAD2DEG { 180.0e0 / DPI };

} // namespace ngpt

#endif
