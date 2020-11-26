///
/// @file meridian_arc_length.hpp
///
/// @brief Algorithms for computation of meridian arc length
///
/// @see https://en.wikipedia.org/wiki/Meridian_arc

#ifndef __MERIDIAN_ARC_LENGTH_HPP__
#define __MERIDIAN_ARC_LENGTH_HPP__

#include "ellipsoid.hpp"

namespace ngpt {

/// @brief core namespace holds the core of ellipsoid-related functions
namespace core {

// @brief Table of coefficients for meridian arc computation
double meridian_arc_length_impl1(double a, double f, double lat) noexcept;
double meridian_arc_length_impl2(double a, double f, double lat) noexcept;
double meridian_arc_length_impl3(double a, double f, double lat) noexcept;
double meridian_arc_length_impl4(double a, double f, double lat) noexcept;

} // namespace core

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

} // namespace ngpt

#endif
