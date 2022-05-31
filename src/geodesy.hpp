/// @file geodesy.hpp
/// @brief A list of frequently used geodetic functions.
/// @see http://www.movable-type.co.uk/scripts/latlong.html

#ifndef __NGPT_GEODESY_HPP__
#define __NGPT_GEODESY_HPP__

#include "geoconst.hpp"
#include "matvec/matvec.hpp"
#include <cassert>
#include <cmath>
#include <stdexcept>
#include <type_traits>

namespace dso {

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
void top2daz(double east, double north, double up, double &distance,
             double &azimouth, double &zenith);
inline
void top2daz(const Vector3 &enu, double &distance, double &azimouth,
             double &zenith) {
  return top2daz(enu.x(), enu.y(), enu.z(), distance, azimouth, zenith);
}

/// @brief Compute distance, azimouth and elevation from topocentric vector
/// @param[in] enu Vector of size >= 3, containing East, North and Up 
///                 coordinates in [m]
/// @param[out] distance  The distance/norm of the topocentric vector [m]
/// @param[out] azimouth  The azimouth between the two points (in the 
///                       topocentric frame) in [rad]. Range [0,2π]
/// @param[out] elevation The elevation between the two points of the vector
///                       in [rad]. Range [0, π]
void top2dae(const Vector3 &enu, double &distance, double &azimouth,
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
void top2dae(const Vector3 &enu, double &distance, double &azimouth,
             double &elevation, Vector3 &dAdr, Vector3 &dEdr);

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
inline
void top2car(const Vector3 &enu, double lat, double lon, Vector3 &dr) noexcept {
  return top2car(enu.x(), enu.y(), enu.z(), lon, lat, dr.x(), dr.y(), dr.z());
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
