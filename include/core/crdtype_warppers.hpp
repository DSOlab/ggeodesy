/** @file
 * Define a number of classes which are actually wrappers around a 3-D vector 
 * to distinguish between Coordinate types.
 */

#ifndef __DSO_COORDINATE_TYPE_WRAPPERS_HPP__
#define __DSO_COORDINATE_TYPE_WRAPPERS_HPP__

#include "eigen3/Eigen/Eigen"

namespace dso {

namespace detail {
  using Vec3d = Eigen::Matrix<double, 3, 1>;
}/* namespace detail */

struct CartesianCrd {
  /* mv = (X,Y,Z) in [m] */
  detail::Vec3d mv;
  double x() const noexcept { return mv(0); }
  double y() const noexcept { return mv(1); }
  double z() const noexcept { return mv(2); }
  double &x() noexcept { return mv(0); }
  double &y() noexcept { return mv(1); }
  double &z() noexcept { return mv(2); }
  CartesianCrd() noexcept : mv{} {};
  explicit CartesianCrd(const detail::Vec3d &vec) noexcept : mv(vec) {};
  CartesianCrd(double x, double y, double z) noexcept {
    mv << x, y, z;
  }
};

struct CartesianCrdView {
  /* mv = (X,Y,Z) in [m] */
  detail::Vec3d& mv;
  explicit CartesianCrdView(detail::Vec3d &v) noexcept : mv(v) {};
  double x() const noexcept { return mv(0); }
  double y() const noexcept { return mv(1); }
  double z() const noexcept { return mv(2); }
  double &x() noexcept { return mv(0); }
  double &y() noexcept { return mv(1); }
  double &z() noexcept { return mv(2); }
};

struct CartesianCrdConstView {
  /* mv = (X,Y,Z) in [m] */
  const detail::Vec3d& mv;
  explicit CartesianCrdConstView(const detail::Vec3d &v) noexcept : mv(v) {};
  double x() const noexcept { return mv(0); }
  double y() const noexcept { return mv(1); }
  double z() const noexcept { return mv(2); }
};

struct GeodeticCrd {
  /* mv = (φ,λ,h) with φ geodetic lattide and h ellipsoidal height. Units
   * in ([rad], [rad], [m]) and ranges:
   * -π/2 <= φ < π/2
   * -π <= λ < π
   * h Real number
   */
  detail::Vec3d mv;
  double lat() const noexcept { return mv(0); }
  double lon() const noexcept { return mv(1); }
  double hgt() const noexcept { return mv(2); }
  double &lat() noexcept { return mv(0); }
  double &lon() noexcept { return mv(1); }
  double &hgt() noexcept { return mv(2); }
};

struct GeodeticCrdView {
  /* mv = (φ,λ,h) with φ geodetic lattide and h ellipsoidal height. Units
   * in ([rad], [rad], [m]) and ranges:
   * -π/2 <= φ < π/2
   * -π <= λ < π
   * h Real number
   */
  detail::Vec3d& mv;
  explicit GeodeticCrdView(detail::Vec3d &v) noexcept : mv(v) {};
  double lat() const noexcept { return mv(0); }
  double lon() const noexcept { return mv(1); }
  double hgt() const noexcept { return mv(2); }
  double &lat() noexcept { return mv(0); }
  double &lon() noexcept { return mv(1); }
  double &hgt() noexcept { return mv(2); }
};

struct GeodeticCrdConstView {
  /* mv = (φ,λ,h) with φ geodetic lattide and h ellipsoidal height. Units
   * in ([rad], [rad], [m]) and ranges:
   * -π/2 <= φ < π/2
   * -π <= λ < π
   * h Real number
   */
  const detail::Vec3d& mv;
  explicit GeodeticCrdConstView(const detail::Vec3d &v) noexcept : mv(v) {};
  double lat() const noexcept { return mv(0); }
  double lon() const noexcept { return mv(1); }
  double hgt() const noexcept { return mv(2); }
};

struct SphericalCrd {
  /* mv = (r,φ,λ) with φ geocentric lattide and r radius. Units
   * in ([rad], [rad], [m]) and ranges:
   * -π/2 <= φ < π/2
   * -π <= λ < π
   * r >= 0
   */
  detail::Vec3d mv;
  double r() const noexcept { return mv(0); }
  double lat() const noexcept { return mv(1); }
  double lon() const noexcept { return mv(2); }
  double &r() noexcept { return mv(0); }
  double &lat() noexcept { return mv(1); }
  double &lon() noexcept { return mv(2); }
};

struct SphericalCrdView {
  /* mv = (r,φ,λ) with φ geocentric lattide and r radius. Units
   * in ([rad], [rad], [m]) and ranges:
   * -π/2 <= φ < π/2
   * -π <= λ < π
   * r >= 0
   */
  detail::Vec3d& mv;
  explicit SphericalCrdView(detail::Vec3d &v) noexcept : mv(v) {};
  double r() const noexcept { return mv(0); }
  double lat() const noexcept { return mv(1); }
  double lon() const noexcept { return mv(2); }
  double &r() noexcept { return mv(0); }
  double &lat() noexcept { return mv(1); }
  double &lon() noexcept { return mv(2); }
};

struct SphericalCrdConstView {
  /* mv = (r,φ,λ) with φ geocentric lattide and r radius. Units
   * in ([rad], [rad], [m]) and ranges:
   * -π/2 <= φ < π/2
   * -π <= λ < π
   * r >= 0
   */
  const detail::Vec3d& mv;
  explicit SphericalCrdConstView(const detail::Vec3d &v) noexcept : mv(v) {};
  double r() const noexcept { return mv(0); }
  double lat() const noexcept { return mv(1); }
  double lon() const noexcept { return mv(2); }
};

template<typename T> struct CoordinateTypeTraits {};

template<>struct CoordinateTypeTraits<CartesianCrd> {
  static constexpr const int isCartesian = true;
  static constexpr const int isConst = false;
};
template<>struct CoordinateTypeTraits<CartesianCrdView> {
  static constexpr const int isCartesian = true;
  static constexpr const int isConst = false;
};
template<>struct CoordinateTypeTraits<CartesianCrdConstView> {
  static constexpr const int isCartesian = true;
  static constexpr const int isConst = true;
};
template<>struct CoordinateTypeTraits<GeodeticCrd> {
  static constexpr const int isGeodetic = true;
  static constexpr const int isConst = false;
};
template<>struct CoordinateTypeTraits<GeodeticCrdView> {
  static constexpr const int isGeodetic = true;
  static constexpr const int isConst = false;
};
template<>struct CoordinateTypeTraits<GeodeticCrdConstView> {
  static constexpr const int isGeodetic = true;
  static constexpr const int isConst = true;
};
template<>struct CoordinateTypeTraits<SphericalCrd> {
  static constexpr const int isSpherical = true;
  static constexpr const int isConst = false;
};
template<>struct CoordinateTypeTraits<SphericalCrdView> {
  static constexpr const int isSpherical = true;
  static constexpr const int isConst = false;
};
template<> struct CoordinateTypeTraits<SphericalCrdConstView> {
  static constexpr const int isSpherical = true;
  static constexpr const int isConst = true;
};

} /* namespace dso */

#endif
