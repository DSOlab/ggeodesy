///
/// @file car2top.hpp
///
/// @brief Transformation of cartesian to topocentric vector.
///
/// The functions here have (in general) two forms, depending on how the user
/// passes in the cartesian geocentric vector \f$\Delta\f$x; that is, users can
/// pass in two points (named i and j with their x, y and z components), or pass
/// in the vector directly (as dx, dy and dz). Note that in this last case, we
/// also need to know the coordinates of the first/starting point (i), cause
/// this is the central point for the topocentric system.
///

#ifndef _CARTESIAN_TO_TOPOCENTRIC_
#define _CARTESIAN_TO_TOPOCENTRIC_

#include "car2ell.hpp"
#include <cmath>

namespace ngpt {

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
  double cosf{std::cos(phi_i)};
  double cosl{std::cos(lambda_i)};
  double sinf{std::sin(phi_i)};
  double sinl{std::sin(lambda_i)};

  // Topocentric vector.
  north = -sinf * cosl * dx - sinf * sinl * dy + cosf * dz;
  east = -sinl * dx + cosl * dy;
  up = cosf * cosl * dx + cosf * sinl * dy + sinf * dz;

  // Finished.
  return;
}

} // namespace core

/// @brief Cartesian to topocentric (vector).
///
/// @tparam     E      The reference ellipsoid (i.e. one of ngpt::ellipsoid).
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
/// @see ngpt::core::dcar2top
template <ellipsoid E>
void car2top(double xi, double yi, double zi, double xj, double yj, double zj,
             double &north, double &east, double &up) noexcept {
  constexpr double semi_major{ellipsoid_traits<E>::a};
  constexpr double flattening{ellipsoid_traits<E>::f};

  // Catresian vector.
  double dx{xj - xi};
  double dy{yj - yi};
  double dz{zj - zi};

  // transform to topocentric
  core::dcar2top(xi, yi, zi, dx, dy, dz, semi_major, flattening, north, east,
                 up);

  // Finished.
  return;
}

/// @brief Cartesian to topocentric (vector).
///
/// @tparam     E      The reference ellipsoid (i.e. one of ngpt::ellipsoid).
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
/// @see ngpt::core::dcar2top
template <ellipsoid E>
void dcar2top(double xi, double yi, double zi, double dx, double dy, double dz,
              double &north, double &east, double &up) noexcept {
  constexpr double semi_major{ellipsoid_traits<E>::a};
  constexpr double flattening{ellipsoid_traits<E>::f};
  core::dcar2top(xi, yi, zi, dx, dy, dz, semi_major, flattening, north, east,
                 up);
}

/// @brief Cartesian to topocentric (vector).
///
/// @param[in]  xi     Cartesian x-component of point i (meters)
/// @param[in]  yi     Cartesian y-component of point i (meters)
/// @param[in]  zi     Cartesian z-component of point i (meters)
/// @param[in]  xj     Cartesian x-component of point j (meters)
/// @param[in]  yj     Cartesian y-component of point j (meters)
/// @param[in]  zj     Cartesian z-component of point j (meters)
/// @param[in]  e      reference ellipsoid (ngpt::Ellipsoid)
/// @param[out] north  Vector north component (meters)
/// @param[out] east   Vector east component (meters)
/// @param[out] up     Vector up component (meters)
///
/// @see ngpt::core::dcar2top
void car2top(double xi, double yi, double zi, double xj, double yj, double zj,
             const Ellipsoid &e, double &north, double &east,
             double &up) noexcept {
  const double semi_major{e.semi_major()};
  const double flattening{e.flattening()};

  // Catresian vector.
  double dx{xj - xi};
  double dy{yj - yi};
  double dz{zj - zi};

  // transform to topocentric
  core::dcar2top(xi, yi, zi, dx, dy, dz, semi_major, flattening, north, east,
                 up);

  // Finished.
  return;
}

/// @brief Cartesian to topocentric (vector).
///
/// @param[in]  xi     Cartesian x-component of point i (meters)
/// @param[in]  yi     Cartesian y-component of point i (meters)
/// @param[in]  zi     Cartesian z-component of point i (meters)
/// @param[in]  dx     x-component of \f$\Delta x\f$ vector (meters)
/// @param[in]  dy     y-component of \f$\Delta y\f$ vector (meters)
/// @param[in]  dz     z-component of \f$\Delta z\f$ vector (meters)
/// @param[in]  e      the reference ellipsoid (ngpt::Ellipsoid)
/// @param[out] north  Vector north component (meters)
/// @param[out] east   Vector east component (meters)
/// @param[out] up     Vector up component (meters)
///
/// @see ngpt::core::dcar2top
void dcar2top(double xi, double yi, double zi, double dx, double dy, double dz,
              const Ellipsoid &e, double &north, double &east,
              double &up) noexcept {
  const double semi_major{e.semi_major()};
  const double flattening{e.flattening()};
  core::dcar2top(xi, yi, zi, dx, dy, dz, semi_major, flattening, north, east,
                 up);
}

} // namespace ngpt

#endif
