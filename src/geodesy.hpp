/** @file
 * A list of frequently used geodetic functions.
 */

#ifndef __DSO_BASIC_GEODESY_HPP__
#define __DSO_BASIC_GEODESY_HPP__

#include "crd_transformations.hpp"

namespace dso {
struct CartesianCrd {
  /* mv = (X,Y,Z) in [m] */
  Eigen::Matrix<double, 3, 1> mv;
  double x() const noexcept { return mv(0); }
  double y() const noexcept { return mv(1); }
  double z() const noexcept { return mv(2); }
  double &x() noexcept { return mv(0); }
  double &y() noexcept { return mv(1); }
  double &z() noexcept { return mv(2); }
};

struct GeodeticCrd {
  /* mv = (φ,λ,h) with φ geodetic lattide and h ellipsoidal height. Units
   * in ([rad], [rad], [m]) and ranges:
   * -π/2 <= φ < π/2
   * -π <= λ < π
   * h Real number
   */
  Eigen::Matrix<double, 3, 1> mv;
  double lat() const noexcept { return mv(0); }
  double lon() const noexcept { return mv(1); }
  double hgt() const noexcept { return mv(2); }
  double &lat() noexcept { return mv(0); }
  double &lon() noexcept { return mv(1); }
  double &hgt() noexcept { return mv(2); }
};

struct SphericalCrd {
  /* mv = (r,φ,λ) with φ geocentric lattide and r radius. Units
   * in ([rad], [rad], [m]) and ranges:
   * -π/2 <= φ < π/2
   * -π <= λ < π
   * r >= 0
   */
  Eigen::Matrix<double, 3, 1> mv;
  double r() const noexcept { return mv(0); }
  double lat() const noexcept { return mv(1); }
  double lon() const noexcept { return mv(2); }
  double &r() noexcept { return mv(0); }
  double &lat() noexcept { return mv(1); }
  double &lon() noexcept { return mv(2); }
};

inline SphericalCrd cartesian2spherical(const CartesianCrd &v) noexcept {
  SphericalCrd s;
  cartesian2spherical(v.x(), v.y(), v.z(), s.r(), s.lat(), s.lon());
  return s;
}

inline CartesianCrd spherical2cartesian(const SphericalCrd &v) noexcept {
  CartesianCrd s;
  spherical2cartesian(v.r(), v.lat(), v.lon(), s.x(), s.y(), s.z());
  return s;
}

template <ellipsoid E>
CartesianCrd geodetic2cartesian(const GeodeticCrd &v) noexcept {
  CartesianCrd s;
  geodetic2cartesian<E>(v.lat(), v.lon(), v.hgt(), s.x(), s.y(), s.z());
  return s;
}

template <ellipsoid E>
GeodeticCrd cartesian2geodetic(const CartesianCrd &v) noexcept {
  GeodeticCrd s;
  cartesian2geodetic<E>(v.x(), v.y(), v.z(), s.lat(), s.lon(), s.hgt());
  return s;
}
} /* namespace dso */

#endif
