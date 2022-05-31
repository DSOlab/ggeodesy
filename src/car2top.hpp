/// @file car2top.hpp
/// @brief Transformation of cartesian to topocentric vector.
/// The functions here have (in general) two forms, depending on how the user
/// passes in the cartesian geocentric vector \f$\Delta\f$x; that is, users can
/// pass in two points (named i and j with their x, y and z components), or pass
/// in the vector directly (as dx, dy and dz). Note that in this last case, we
/// also need to know the coordinates of the first/starting point (i), cause
/// this is the central point for the topocentric system.

#ifndef __CARTESIAN_TO_TOPOCENTRIC__
#define __CARTESIAN_TO_TOPOCENTRIC__

#include "car2ell.hpp"
#include "matvec/matvec.hpp"
#include <cmath>

namespace dso {

namespace core {

/// @brief Cartesian to topocentric (vector).
///
/// Transform a vector expressed in cartesian, geocentric coordinates to the
/// topocentric, local system around point i. This function depends on the
/// reference ellipsoid. All units in meters/radians. The function will
/// transform the (geocentric) vector \f$\vec{\Delta X}\f$ to the local
/// topocentric reference frame around point \f$\vec{X}_i\f$.
///
/// @param[in]  xi     Cartesian x-component of point i (meters)
/// @param[in]  yi     Cartesian y-component of point i (meters)
/// @param[in]  zi     Cartesian z-component of point i (meters)
/// @param[in]  dx     x-component of \f$\Delta x\f$ vector (meters)
/// @param[in]  dy     y-component of \f$\Delta y\f$ vector (meters)
/// @param[in]  dz     z-component of \f$\Delta z\f$ vector (meters)
/// @param[in]  semi_major  The semi-major axis of the ref. ellipsoid
/// @param[in]  flattening  The flattening of the ref. ellipsoid
/// @param[out] north  Vector north component (meters)
/// @param[out] east   Vector east component (meters)
/// @param[out] up     Vector up component (meters)
/// @throw             Does not throw.
///
/// @note The ellispoid is needed to transform the cartesian coordinates of
///       the (reference) point i to ellispoidal coordinates.
///
/// @see "Physical Geodesy", pg. 209
///
void dcar2top(double xi, double yi, double zi, double dx, double dy, double dz,
              double semi_major, double flattening, double &north, double &east,
              double &up) noexcept {

  // Ellipsoidal coordinates of reference point.
  double phi_i, lambda_i, h_i;

  // Cartesian to ellipsoidal for reference point.
  core::car2ell(xi, yi, zi, semi_major, flattening, phi_i, lambda_i, h_i);

  // Trigonometric numbers.
  const double cosf{std::cos(phi_i)};
  const double cosl{std::cos(lambda_i)};
  const double sinf{std::sin(phi_i)};
  const double sinl{std::sin(lambda_i)};

  // Topocentric vector.
  north = -sinf * cosl * dx - sinf * sinl * dy + cosf * dz;
  east = -sinl * dx + cosl * dy;
  up = cosf * cosl * dx + cosf * sinl * dy + sinf * dz;

  // Finished.
  return;
}

Vector3 dcar2top(const Vector3 &r, const Vector3 &dr, double semi_major,
                 double flattening) noexcept {

  // Cartesian to ellipsoidal for reference point.
  const Vector3 lfh = core::car2ell(r, semi_major, flattening);

  // Trigonometric numbers.
  const double cosf{std::cos(lfh.y())};
  const double sinf{std::sin(lfh.y())};
  const double sinl{std::sin(lfh.x())};
  const double cosl{std::cos(lfh.x())};

  // Topocentric vector.
  const double north =
      -sinf * cosl * dr.x() - sinf * sinl * dr.y() + cosf * dr.z();
  const double east = -sinl * dr.x() + cosl * dr.y();
  const double up = cosf * cosl * dr.x() + cosf * sinl * dr.y() + sinf * dr.z();

  // Finished.
  return Vector3({east, north, up});
}

} // namespace core

/// @brief Cartesian to topocentric (vector).
///
/// @tparam     E      The reference ellipsoid (i.e. one of dso::ellipsoid).
/// @param[in]  xi     Cartesian x-component of point i (meters)
/// @param[in]  yi     Cartesian y-component of point i (meters)
/// @param[in]  zi     Cartesian z-component of point i (meters)
/// @param[in]  xj     Cartesian x-component of point j (meters)
/// @param[in]  yj     Cartesian y-component of point j (meters)
/// @param[in]  zj     Cartesian z-component of point j (meters)
/// @param[out] north  Vector north component (meters)
/// @param[out] east   Vector east component (meters)
/// @param[out] up     Vector up component (meters)
///
/// @see dso::core::dcar2top
template <ellipsoid E>
void car2top(double xi, double yi, double zi, double xj, double yj, double zj,
             double &north, double &east, double &up) noexcept {
  constexpr double semi_major{ellipsoid_traits<E>::a};
  constexpr double flattening{ellipsoid_traits<E>::f};

  // Catresian vector.
  const double dx{xj - xi};
  const double dy{yj - yi};
  const double dz{zj - zi};

  // transform to topocentric
  core::dcar2top(xi, yi, zi, dx, dy, dz, semi_major, flattening, north, east,
                 up);

  // Finished.
  return;
}
template <ellipsoid E>
Vector3 car2top(const Vector3 &xyz_i, const Vector3 &xyz_j) noexcept {
  constexpr double semi_major{ellipsoid_traits<E>::a};
  constexpr double flattening{ellipsoid_traits<E>::f};

  // Catresian vector.
  const Vector3 dr = xyz_j - xyz_i;

  // transform to topocentric
  return core::dcar2top(xyz_i, dr, semi_major, flattening);
}

/// @brief Cartesian to topocentric (vector).
///
/// @tparam     E      The reference ellipsoid (i.e. one of dso::ellipsoid).
/// @param[in]  xi     Cartesian x-component of point i (meters)
/// @param[in]  yi     Cartesian y-component of point i (meters)
/// @param[in]  zi     Cartesian z-component of point i (meters)
/// @param[in]  dx     x-component of \f$\Delta x\f$ vector (meters)
/// @param[in]  dy     y-component of \f$\Delta y\f$ vector (meters)
/// @param[in]  dz     z-component of \f$\Delta z\f$ vector (meters)
/// @param[out] north  Vector north component (meters)
/// @param[out] east   Vector east component (meters)
/// @param[out] up     Vector up component (meters)
/// @throw             Does not throw.
///
/// @see dso::core::dcar2top
template <ellipsoid E>
void dcar2top(double xi, double yi, double zi, double dx, double dy, double dz,
              double &north, double &east, double &up) noexcept {
  constexpr double semi_major{ellipsoid_traits<E>::a};
  constexpr double flattening{ellipsoid_traits<E>::f};
  core::dcar2top(xi, yi, zi, dx, dy, dz, semi_major, flattening, north, east,
                 up);
}
template <ellipsoid E>
Vector3 dcar2top(const Vector3 &xyz_i, const Vector3 &dr) noexcept {
  constexpr double semi_major{ellipsoid_traits<E>::a};
  constexpr double flattening{ellipsoid_traits<E>::f};
  return core::dcar2top(xyz_i, dr, semi_major, flattening);
}

/// @brief Cartesian to topocentric (vector).
///
/// @param[in]  xi     Cartesian x-component of point i (meters)
/// @param[in]  yi     Cartesian y-component of point i (meters)
/// @param[in]  zi     Cartesian z-component of point i (meters)
/// @param[in]  xj     Cartesian x-component of point j (meters)
/// @param[in]  yj     Cartesian y-component of point j (meters)
/// @param[in]  zj     Cartesian z-component of point j (meters)
/// @param[in]  e      reference ellipsoid (dso::Ellipsoid)
/// @param[out] north  Vector north component (meters)
/// @param[out] east   Vector east component (meters)
/// @param[out] up     Vector up component (meters)
///
/// @see dso::core::dcar2top
void car2top(double xi, double yi, double zi, double xj, double yj, double zj,
             const Ellipsoid &e, double &north, double &east,
             double &up) noexcept {
  const double semi_major{e.semi_major()};
  const double flattening{e.flattening()};

  // Catresian vector.
  const double dx{xj - xi};
  const double dy{yj - yi};
  const double dz{zj - zi};

  // transform to topocentric
  core::dcar2top(xi, yi, zi, dx, dy, dz, semi_major, flattening, north, east,
                 up);

  // Finished.
  return;
}
Vector3 car2top(const Vector3 &xyz_i, const Vector3 &xyz_j,
                const Ellipsoid &e) noexcept {
  const double semi_major{e.semi_major()};
  const double flattening{e.flattening()};

  // transform to topocentric
  return core::dcar2top(xyz_i, (xyz_j - xyz_i), semi_major, flattening);
}

/// @brief Cartesian to topocentric (vector).
///
/// @param[in]  xi     Cartesian x-component of point i (meters)
/// @param[in]  yi     Cartesian y-component of point i (meters)
/// @param[in]  zi     Cartesian z-component of point i (meters)
/// @param[in]  dx     x-component of \f$\Delta x\f$ vector (meters)
/// @param[in]  dy     y-component of \f$\Delta y\f$ vector (meters)
/// @param[in]  dz     z-component of \f$\Delta z\f$ vector (meters)
/// @param[in]  e      the reference ellipsoid (dso::Ellipsoid)
/// @param[out] north  Vector north component (meters)
/// @param[out] east   Vector east component (meters)
/// @param[out] up     Vector up component (meters)
///
/// @see dso::core::dcar2top
void dcar2top(double xi, double yi, double zi, double dx, double dy, double dz,
              const Ellipsoid &e, double &north, double &east,
              double &up) noexcept {
  const double semi_major{e.semi_major()};
  const double flattening{e.flattening()};
  core::dcar2top(xi, yi, zi, dx, dy, dz, semi_major, flattening, north, east,
                 up);
}
Vector3 dcar2top(const Vector3 &xyz_i, const Vector3 &dr,
                 const Ellipsoid &e) noexcept {
  const double semi_major{e.semi_major()};
  const double flattening{e.flattening()};
  return core::dcar2top(xyz_i, dr, semi_major, flattening);
}

} // namespace dso

#endif
