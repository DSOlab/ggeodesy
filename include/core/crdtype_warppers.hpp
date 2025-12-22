/** @file
 * Define a number of classes which are actually wrappers around a 3-D vector
 * to distinguish between Coordinate types.
 */

#ifndef __DSO_COORDINATE_TYPE_WRAPPERS_CORE_HPP__
#define __DSO_COORDINATE_TYPE_WRAPPERS_CORE_HPP__

#include "eigen3/Eigen/Eigen"

namespace dso {

namespace detail {
using Vec3d = Eigen::Matrix<double, 3, 1>;
} /* namespace detail */

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
  CartesianCrd(double x, double y, double z) noexcept { mv << x, y, z; }

  const detail::Vec3d &const_ref_vec3d() const noexcept {return mv;}
  detail::Vec3d copy_vec3d() const noexcept {return mv;}
};

struct CartesianCrdView {
  /* mv = (X,Y,Z) in [m] */
  detail::Vec3d &mv;

  explicit CartesianCrdView(detail::Vec3d &v) noexcept : mv(v) {};
  explicit CartesianCrdView(CartesianCrd &v) noexcept : mv(v.mv) {};

  double x() const noexcept { return mv(0); }
  double y() const noexcept { return mv(1); }
  double z() const noexcept { return mv(2); }
  double &x() noexcept { return mv(0); }
  double &y() noexcept { return mv(1); }
  double &z() noexcept { return mv(2); }
  
  const detail::Vec3d &const_ref_vec3d() const noexcept {return mv;}
  detail::Vec3d copy_vec3d() const noexcept {return detail::Vec3d(mv);}
};

struct CartesianCrdConstView {
  /* mv = (X,Y,Z) in [m] */
  const detail::Vec3d &mv;

  explicit CartesianCrdConstView(const detail::Vec3d &v) noexcept : mv(v) {};
  explicit CartesianCrdConstView(const CartesianCrd &v) noexcept : mv(v.mv) {};
  explicit CartesianCrdConstView(const CartesianCrdView &v) noexcept
      : mv(v.mv) {};

  double x() const noexcept { return mv(0); }
  double y() const noexcept { return mv(1); }
  double z() const noexcept { return mv(2); }
  
  const detail::Vec3d &const_ref_vec3d() const noexcept {return mv;}
  detail::Vec3d copy_vec3d() const noexcept {return detail::Vec3d(mv);}
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
  detail::Vec3d &mv;

  explicit GeodeticCrdView(detail::Vec3d &v) noexcept : mv(v) {};
  explicit GeodeticCrdView(GeodeticCrd &v) noexcept : mv(v.mv) {};

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
  const detail::Vec3d &mv;

  explicit GeodeticCrdConstView(const detail::Vec3d &v) noexcept : mv(v) {};
  explicit GeodeticCrdConstView(const GeodeticCrd &v) noexcept : mv(v.mv) {};
  explicit GeodeticCrdConstView(const GeodeticCrdView &v) noexcept
      : mv(v.mv) {};

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
  detail::Vec3d &mv;

  explicit SphericalCrdView(detail::Vec3d &v) noexcept : mv(v) {};
  explicit SphericalCrdView(SphericalCrd &v) noexcept : mv(v.mv) {};

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
  const detail::Vec3d &mv;

  explicit SphericalCrdConstView(const detail::Vec3d &v) noexcept : mv(v) {};
  explicit SphericalCrdConstView(const SphericalCrd &v) noexcept : mv(v.mv) {};
  explicit SphericalCrdConstView(const SphericalCrdView &v) noexcept
      : mv(v.mv) {};

  double r() const noexcept { return mv(0); }
  double lat() const noexcept { return mv(1); }
  double lon() const noexcept { return mv(2); }
};

/** Generic traits for all Coordinate triplets (to be specialized...) */
template <typename T> struct CoordinateTypeTraits {};

template <> struct CoordinateTypeTraits<CartesianCrd> {
  static constexpr const int isCartesian = true;
  static constexpr const int isConst = false;
#if !defined(__cplusplus) || __cplusplus < 202002L
    /* We're in C++17 or earlier, no concepts */
    static constexpr const int isGeodetic  = false;
    static constexpr const int isSpherical = false;
#endif
};

template <> struct CoordinateTypeTraits<CartesianCrdView> {
  static constexpr const int isCartesian = true;
  static constexpr const int isConst = false;
#if !defined(__cplusplus) || __cplusplus < 202002L
    static constexpr const int isGeodetic  = false;
    static constexpr const int isSpherical = false;
#endif
};

template <> struct CoordinateTypeTraits<CartesianCrdConstView> {
  static constexpr const int isCartesian = true;
  static constexpr const int isConst = true;
#if !defined(__cplusplus) || __cplusplus < 202002L
    static constexpr const int isGeodetic  = false;
    static constexpr const int isSpherical = false;
#endif
};

template <> struct CoordinateTypeTraits<GeodeticCrd> {
  static constexpr const int isGeodetic = true;
  static constexpr const int isConst = false;
#if !defined(__cplusplus) || __cplusplus < 202002L
    static constexpr const int isCartesian = false;
    static constexpr const int isSpherical = false;
#endif
};

template <> struct CoordinateTypeTraits<GeodeticCrdView> {
  static constexpr const int isGeodetic = true;
  static constexpr const int isConst = false;
#if !defined(__cplusplus) || __cplusplus < 202002L
    static constexpr const int isCartesian = false;
    static constexpr const int isSpherical = false;
#endif
};

template <> struct CoordinateTypeTraits<GeodeticCrdConstView> {
  static constexpr const int isGeodetic = true;
  static constexpr const int isConst = true;
#if !defined(__cplusplus) || __cplusplus < 202002L
    static constexpr const int isCartesian = false;
    static constexpr const int isSpherical = false;
#endif
};

template <> struct CoordinateTypeTraits<SphericalCrd> {
  static constexpr const int isSpherical = true;
  static constexpr const int isConst = false;
#if !defined(__cplusplus) || __cplusplus < 202002L
    static constexpr const int isCartesian = false;
    static constexpr const int isGeodetic  = false;
#endif
};

template <> struct CoordinateTypeTraits<SphericalCrdView> {
  static constexpr const int isSpherical = true;
  static constexpr const int isConst = false;
#if !defined(__cplusplus) || __cplusplus < 202002L
    static constexpr const int isCartesian = false;
    static constexpr const int isGeodetic  = false;
#endif
};

template <> struct CoordinateTypeTraits<SphericalCrdConstView> {
  static constexpr const int isSpherical = true;
  static constexpr const int isConst = true;
#if !defined(__cplusplus) || __cplusplus < 202002L
    /* We're in C++17 or earlier, no concepts */
    static constexpr const int isCartesian = false;
    static constexpr const int isGeodetic  = false;
#endif
};


/* Concepts for C++20 upwards */
#if __cplusplus >= 202002L

/* The following two templates, enabling matching types T that have a member 
 * function Traits<T>::isCartesian
 */
template<typename, typename = void> struct has_is_cartesian : std::false_type {};
template<typename T>
struct has_is_cartesian<T, std::void_t<decltype(Traits<T>::isCartesian)>> : std::true_type {};

/* Now we define a concept: Any type T that :
 * 1. has a Traits<T>::isCartesian, and
 * 2. the value of Traits<T>::isCartesian is true
 */
template<typename T> concept IsCartesian = 
  has_is_cartesian<T>::value && Traits<T>::isCartesian;

/* The following two templates, enabling matching types T that have a member 
 * function Traits<T>::isSpherical
 */
template<typename, typename = void> struct has_is_spherical : std::false_type {};
template<typename T>
struct has_is_spherical<T, std::void_t<decltype(Traits<T>::isSpherical)>> : std::true_type {};

/* Now we define a concept: Any type T that :
 * 1. has a Traits<T>::isSpherical, and
 * 2. the value of Traits<T>::isSpherical is true
 */
template<typename T> concept IsSpherical = 
  has_is_cartesian<T>::value && Traits<T>::isSpherical;

/* The following two templates, enabling matching types T that have a member 
 * function Traits<T>::isGeodetic
 */
template<typename, typename = void> struct has_is_geodetic : std::false_type {};
template<typename T>
struct has_is_geodetic<T, std::void_t<decltype(Traits<T>::isGeodetic)>> : std::true_type {};

/* Now we define a concept: Any type T that :
 * 1. has a Traits<T>::isGeodetic, and
 * 2. the value of Traits<T>::isGeodetic is true
 */
template<typename T> concept IsGeodetic = 
  has_is_cartesian<T>::value && Traits<T>::isGeodetic;

#endif

} /* namespace dso */

#endif
