/// @file geodesy.hpp
/// @brief A list of frequently used geodetic functions.
/// @see http://www.movable-type.co.uk/scripts/latlong.html

#ifndef __NGPT_GEODESY_HPP__
#define __NGPT_GEODESY_HPP__

#include "eigen3/Eigen/Eigen"
#include "ellipsoid.hpp"
#include "geoconst.hpp"
#include <cassert>
#include <cmath>
#include <stdexcept>
#include <type_traits>

namespace dso {

/// @brief Given the geodetic coordinates of a reference point, return the
///        matrix that turns any vector from the ference point to point P to
///        topocentric coordinates (e,n,u)
Eigen::Matrix<double, 3, 3> topocentric_matrix(double lambda,
                                               double phi) noexcept;

inline Eigen::Matrix<double, 3, 3>
topocentric_matrix(const Eigen::Matrix<double, 3, 1> &lfh) noexcept {
  return topocentric_matrix(lfh(0), lfh(1));
}

/// @brief Ellipsoidal to cartesian coordinates.
///
/// Transform (geocentric) cartesian coordinates (on the ellipsoid) to
/// ellipsoidal coordinates. Units are meters and radians.
///
/// @tparam      E      The reference ellipsoid (i.e. one of dso::ellipsoid).
/// @param[in]   phi    Ellipsoidal latitude (radians)
/// @param[in]   lambda Ellipsoidal longtitude (radians)
/// @param[in]   h      Ellipsoidal height (meters)
/// @param[out]  x      Cartesian x-component (meters)
/// @param[out]  y      Cartesian y-component (meters)
/// @param[out]  z      Cartesian z-component (meters)
/// @throw              Does not throw.
///
template <ellipsoid E>
void ell2car(double lambda, double phi, double h, double &x, double &y,
             double &z) noexcept {
  // Eccentricity squared.
  constexpr double e2{dso::eccentricity_squared<E>()};

  // Radius of curvature in the prime vertical.
  const double N{dso::N<E>(phi)};

  // Trigonometric numbers.
  const double sinf{std::sin(phi)};
  const double cosf{std::cos(phi)};
  const double sinl{std::sin(lambda)};
  const double cosl{std::cos(lambda)};

  // Compute geocentric rectangular coordinates.
  x = (N + h) * cosf * cosl;
  y = (N + h) * cosf * sinl;
  z = ((1e0 - e2) * N + h) * sinf;

  // Finished.
  return;
}

template <ellipsoid E>
Eigen::Matrix<double, 3, 1>
ell2car(const Eigen::Matrix<double, 3, 1> &lfh) noexcept {
  // Eccentricity squared.
  constexpr double e2{dso::eccentricity_squared<E>()};

  // Radius of curvature in the prime vertical.
  const double N{dso::N<E>(lfh(1))};

  // Trigonometric numbers.
  const double sinf{std::sin(lfh(1))};
  const double cosf{std::cos(lfh(1))};
  const double sinl{std::sin(lfh(0))};
  const double cosl{std::cos(lfh(0))};

  // Compute geocentric rectangular coordinates.
  const double x = (N + lfh(2)) * cosf * cosl;
  const double y = (N + lfh(2)) * cosf * sinl;
  const double z = ((1e0 - e2) * N + lfh(2)) * sinf;

// Finished.
  /*const*/ double data[] = {x, y, z};
  return Eigen::Map<Eigen::Matrix<double, 3, 1>>(data, 3);
}

void ell2car(double lambda, double phi, double h, const Ellipsoid &e, double &x,
             double &y, double &z) noexcept;
Eigen::Matrix<double, 3, 1> ell2car(const Eigen::Matrix<double, 3, 1> &lfh,
                                    const Ellipsoid &e) noexcept;

/// @brief Topocentric vector to azimouth, zenith and distance.
///
/// Compute the Distance, Azimouth and Zenith distance given a vector expressed
/// in a local, topocentric system (i.e. given the north, east and up
/// components of the vector).
///
/// @param[in]  north    Vector north component (meters)
/// @param[in]  east     Vector east component (meters)
/// @param[in]  up       Vector up component (meters)
/// @param[out] distance The length of the vector (meters)
/// @param[out] azimouth The azimouth of the vector, radians [0,2*pi).
/// @param[out] zenith   The zenith distance, radians [0,pi)
/// @throw               std::runtime_error if zero division encountered.
///
/// @see "Physical Geodesy", pg. 210
///
double top2daz(double east, double north, double up, double &azimouth,
               double &zenith);

inline double top2daz(const Eigen::Matrix<double, 3, 1> &enu, double &azimouth,
                      double &zenith) {
  return top2daz(enu(0), enu(1), enu(2), azimouth, zenith);
}

/// @brief Compute distance, azimouth and elevation from topocentric vector
/// @param[in] enu Vector of size >= 3, containing East, North and Up
///                 coordinates in [m]
/// @param[out] distance  The distance/norm of the topocentric vector [m]
/// @param[out] azimouth  The azimouth between the two points (in the
///                       topocentric frame) in [rad]. Range [0,2π]
/// @param[out] elevation The elevation between the two points of the vector
///                       in [rad]. Range [0, π]
double top2dae(const Eigen::Matrix<double, 3, 1> &enu, double &azimouth,
               double &elevation);

/// @brief Compute distance, azimouth and elevation from topocentric vector and
///        their partial derivatives w.r.t the topocentic RF
/// @param[in] enu Vector of size >= 3, containing East, North and Up
///                 coordinates in [m]
/// @param[out] distance  The distance/norm of the topocentric vector [m]
/// @param[out] azimouth  The azimouth between the two points (in the
///                       topocentric frame) in [rad]. Range [0,2π]
/// @param[out] elevation The elevation between the two points of the vector
///                       in [rad]. Range [0, π]
/// @param[out] dAdr      Partials of the Azimouth, w.r.t the [e,n,u] (unit)
///                       vectors
/// @param[out] dedr      Partials of the Elevation, w.r.t the [e,n,u] (unit)
///                       vectors
/// @see Satellite Orbits: Models, Methods and Applications, ch 7.4
double top2dae(const Eigen::Matrix<double, 3, 1> &enu, double &azimouth,
               double &elevation, Eigen::Matrix<double, 3, 1> &dAdr,
               Eigen::Matrix<double, 3, 1> &dEdr);

/// @brief Topocentric vector to a geocentric, cartesian vector
///
/// Given a vector expressed in a local, topocentric system (i.e. given the
/// north, east and up components of the vector) around point (lon, lat),
/// transform the vector to a geocentric, cartesian one.
///
/// @param[in]  north    Vector north component (meters)
/// @param[in]  east     Vector east component (meters)
/// @param[in]  up       Vector up component (meters)
/// @param[in]  lat      Latitude of reference point of the vector (radians)
/// @param[in]  lon      Longtitude of reference point of the vector (radians)
/// @param[out] dx       x-component of the cartesian vector (meters)
/// @param[out] dy       y-component of the cartesian vector (meters)
/// @param[out] dz       z-component of the cartesian vector (meters)
///
/// @see "Physical Geodesy", pg. 210
///
void top2car(double east, double north, double up, double lon, double lat,
             double &dx, double &dy, double &dz) noexcept;
inline Eigen::Matrix<double, 3, 1>
top2car(const Eigen::Matrix<double, 3, 1> &enu, double lon,
        double lat) noexcept {
  return topocentric_matrix(lon, lat) * enu;
}
inline Eigen::Matrix<double, 3, 1>
top2car(const Eigen::Matrix<double, 3, 1> &enu,
        const Eigen::Matrix<double, 3, 1> &lfh) noexcept {
  return top2car(enu, lfh(0), lfh(1));
}

/// @brief Bearing (i.e. forward azimouth) of great circle between two points
/// on the sphere.
///
/// This formula is for the initial bearing (sometimes referred to as forward
/// azimuth) which if followed in a straight line along a great-circle arc will
/// take you from the start point to the end point.
///
/// @tparam     T     Any floating point type for input and results.
/// @param[in]  lat1  Latitude of starting point in radians.
/// @param[in]  lon1  Longtitude of starting point in radians.
/// @param[in]  lat2  Latitude of ending point in radians.
/// @param[in]  lon2  Longtitude of ending point in radians.
/// @throw            Does not throw
///
/// @bug This give way too big a difference from Vincenty.
template <typename T,
          typename = std::enable_if_t<std::is_floating_point<T>::value>>
T bearing(T lat1, T lon1, T lat2, T lon2) noexcept {
  using std::cos;
  using std::sin;

  const double DeltaLambda = lon2 - lon1;
  const double cosLat2 = cos(lat2);
  const double nom = sin(DeltaLambda) * cosLat2;
  const double denom =
      cos(lat1) * sin(lat2) - sin(lat1) * cosLat2 * cos(DeltaLambda);
  const double angle = std::atan2(nom, denom);
  // normalize to [0-2pi)
  return std::fmod(angle + D2PI, D2PI);
}

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
              double semi_major, double flattening, double &east, double &north,
              double &up) noexcept;
Eigen::Matrix<double, 3, 1> dcar2top(const Eigen::Matrix<double, 3, 1> &r,
                                     const Eigen::Matrix<double, 3, 1> &dr,
                                     double semi_major,
                                     double flattening) noexcept;
void car2ell(double x, double y, double z, double semi_major, double flattening,
             double &lambda, double &phi, double &h) noexcept;

} // namespace core

Eigen::Matrix<double, 3, 1> car2ell(const Eigen::Matrix<double, 3, 1> &xyz,
                                    double semi_major,
                                    double flattening) noexcept;

/// @brief Cartesian to Spherical, given the radius
/// See Physical Geodesy, Section 1.4
/// The returned vector is: [r, θ, λ]
/// with r: radius vector, θ polar distance and λ geocentric longitude
inline Eigen::Matrix<double, 3, 1>
car2sph(const Eigen::Matrix<double, 3, 1> &xyz, double R) noexcept {
  const double r = xyz.norm();
  const double x = xyz(0);
  const double y = xyz(1);
  const double z = xyz(2);
  const double theta = std::atan2(std::sqrt(x * x + y * y), z);
  const double lambda = std::atan2(y, x);
  return Eigen::Matrix<double, 3, 1>(r, theta, lambda);
}

template <ellipsoid E>
Eigen::Matrix<double, 3, 1>
car2ell(const Eigen::Matrix<double, 3, 1> &xyz) noexcept {
  constexpr double semi_major{ellipsoid_traits<E>::a};
  constexpr double flattening{ellipsoid_traits<E>::f};
  return car2ell(xyz, semi_major, flattening);
}

inline Eigen::Matrix<double, 3, 1>
car2ell(const Eigen::Matrix<double, 3, 1> &xyz, const Ellipsoid &e) noexcept {
  double semi_major{e.semi_major()};
  double flattening{e.flattening()};
  return car2ell(xyz, semi_major, flattening);
}

inline Eigen::Matrix<double, 3, 1>
car2ell(const Eigen::Matrix<double, 3, 1> &xyz, ellipsoid e) noexcept {
  return car2ell(xyz, Ellipsoid(e));
}

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
Eigen::Matrix<double, 3, 1>
car2top(const Eigen::Matrix<double, 3, 1> &xyz_i,
        const Eigen::Matrix<double, 3, 1> &xyz_j) noexcept {
  constexpr double semi_major{ellipsoid_traits<E>::a};
  constexpr double flattening{ellipsoid_traits<E>::f};

  // Catresian vector.
  const Eigen::Matrix<double, 3, 1> dr = xyz_j - xyz_i;

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
Eigen::Matrix<double, 3, 1>
dcar2top(const Eigen::Matrix<double, 3, 1> &xyz_i,
         const Eigen::Matrix<double, 3, 1> &dr) noexcept {
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
inline Eigen::Matrix<double, 3, 1>
car2top(const Eigen::Matrix<double, 3, 1> &xyz_i,
        const Eigen::Matrix<double, 3, 1> &xyz_j, const Ellipsoid &e) noexcept {
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
inline Eigen::Matrix<double, 3, 1>
dcar2top(const Eigen::Matrix<double, 3, 1> &xyz_i,
         const Eigen::Matrix<double, 3, 1> &dr, const Ellipsoid &e) noexcept {
  const double semi_major{e.semi_major()};
  const double flattening{e.flattening()};
  return core::dcar2top(xyz_i, dr, semi_major, flattening);
}

/// @brief Transformation parameters for PZ-90 to WGS84 reference frames
///
/// This struct holds transformation parameters for a 7-parameter conversion
/// between PZ90 and WGS84 reference frames. Various parameters sets can be
/// used, and most are discussed here:
/// ITRS, PZ-90 and WGS 84: current realizationsand the related transformation
/// parameters, C. Boucher, Z. Altamimi, Journal of Geodesy, November 2001
///
/// The general transformation of the Cartesian coordinates (X) of any point
/// close to the Earth from any one TRS to any other one will be givenby a
/// tri-dimensional similarity transformation (T is a translation vector,
/// \f$\lambda\f$ a scale factor and R a rotation matrix) as follows: \f$X_trs1
/// = T + \lambda * R * X_trs2\f$ The parameters for transforming an X system
/// into an XS system are denoted T1,T2,T3,D,R1,R2, and R3:
///
/// | Xs |   | X |   | T1 |   | D   -R3  R2 | | X |
/// | Ys | = | Y | + | T2 | + | R3   D  -R1 |*| Y |
/// | Zs |   | Z |   | T3 |   |-R2   R1  D  | | Z |
///
constexpr struct {
  double tx, ty, tz, // meters
      r1, r2, r3,    // mas
      d;             // ppb
} pz2wgs_parameters[] = {{7e-2, 0e0, -77e-2, -19e0, -4e0, 353e0, -3e0},
                         {0e0, 250e-2, 0e0, 0e0, 0e0, 392e0, 0e0},
                         {0e0, 0e0, 0e0, 0e0, 0e0, 330e0, 0e0},
                         {-47e-2, -51e-2, -200e-2, 2e0, 1e0, 356e0, 22e0},
                         {-110e-2, -30e-2, -90e-2, 0e0, 0e0, 169e0, -120e0},
                         {0e0, 0e0, -110e-2, -16e0, -4e0, 357e0, 9e0},
                         {-3e-2, -2e-2, -45e-2, -37e0, 10e0, 350e0, 13e0},
                         {30e-2, -10e-2, -90e-2, -3e0, -13e0, 355e0, 0e0},
                         {24e-2, -15e-2, -77e-2, -3e0, -19e0, 353e0, -31e0}};

/// @brief Transform WGS84 to PZ90 coordinates
/// @param[in] xwgs A set of 3-dimensional coordinates in WGS84 rf in meters.
///                 The size of this array should be 3*pts, aka more than one
///                 points can be passed in the following order:
///                 [x1,y1,z1, z2,y2,z2, ... ,xpts,ypts,zpts]
/// @param[out] xpz The resulting input coordinates in the PZ90 rf. The size of
///                 this array should be the same as xwgs (aka 3*pts) and the
///                 order will be the same as in the input array. All units are
///                 meters
/// @param[in] pts  Number of points in the input array
/// @param[in] selection The selection of transformation parameters; the
///                 parameters chosen by this function are :
///                 wgs2pz_parameters[selection]
void pz90_to_wgs84(const double *xwgs, double *xpz, int pts = 1,
                   int selection = 0);

} // namespace dso

#endif
