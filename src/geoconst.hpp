/// @file geoconst.hpp
/// @brief Fundamental constants frequently used within geodetic calculations.

#ifndef __NGPT_GEOCONST_HPP__
#define __NGPT_GEOCONST_HPP__

#include <cmath>

namespace dso {
/// The value of \f$\pi\f$.
#if defined(__GNUC__) && !defined(__llvm__)
constexpr double DPI{std::atan(1e0) * 4e0};
#else
#define _USE_MATH_DEFINES
constexpr double DPI{M_PI};
#endif

/// The value of \f$2*\pi\f$.
constexpr double D2PI{2e0 * DPI};

/// Degrees to Radians coefficient.
constexpr double DEG2RAD{DPI / 180e0};

/// Radians to Degrees coefficient.
constexpr double RAD2DEG{180e0 / DPI};

/// mas to radians factor aka \f$\theta_rad = \theta_mas * MAS2RAD\f$
constexpr double MAS2RAD{4.847309743e-9};

} // namespace dso

#endif
