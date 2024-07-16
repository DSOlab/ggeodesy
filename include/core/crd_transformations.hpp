/** @file
 * A list of frequently used coordinate transformation.
 */

#ifndef __DSO_COORDINATE_TRANSFORMATIONS_HPP__
#define __DSO_COORDINATE_TRANSFORMATIONS_HPP__

#include "geoconst.hpp"
#include "eigen3/Eigen/Eigen"

namespace dso {

/** @brief Cartesian to spherical (geographic) coordinates.
 *
 * @param[in]  x    Cartesian x-component [m]
 * @param[in]  y    Cartesian y-component [m]
 * @param[in]  z    Cartesian z-component [m]
 * @param[out] r    Radius (center of spheriod to point) [m]
 * @param[out] lat  Geocentric latitude [rad] (i.e. 90[deg]-polar distance)
 *                  in range [-π/2, π/2)
 * @param[out] lon  Longitude [rad] in range [-π, π)
 */
void cartesian2spherical(double x, double y, double z, double &r, double &glat,
                         double &lon) noexcept;

/** @brief Spherical (geographic) to Cartesian coordinates.
 *
 * @param[in] r    Radius (center of spheriod to point) [m]
 * @param[in] lat  Geocentric latitude [rad] (i.e. 90[deg]-polar distance)
 *                 in range [-π/2, π/2)
 * @param[in] lon  Longitude [rad] in range [-π, π)
 * @param[out]  x  Cartesian x-component [m]
 * @param[out]  y  Cartesian y-component [m]
 * @param[out]  z  Cartesian z-component [m]
 */
void spherical2cartesian(double r, double glat, double lon, double &x,
                         double &y, double &z) noexcept;

/** Given a point on the ellipsoid/spheroid with (φ,λ), compute topocentric 
 *  rotation matrix.
 *
 * This function computes the unit topocentric vectors (often also called 
 * Local-Vertical Local Tangent) given a reference point on the 
 * ellipsoid or spheroid anc concatenated them in a rotation matrix R. 
 * If (φ,λ) are ellipsoidal coordinates, the  vector u is orthogonal to the
 * tangent plane to the ellipsoid, which is defined by (e,n). If (φ,λ) are 
 * taken as the spherical longitude and latitude, thence, the vector u is in 
 * the radial direction and defines the tangent plane to the sphere.
 *
 * The resulting matrix from this function, R, can be used in the sense:
 * |δX|        |e|
 * |δY|  = R * |n|
 * |δZ|        |u|
 * where (X, Y, Z) are the coordinates in the Cartesian system and (e,n,u) are 
 * coordinates in the East, North and Up directions; the provided (φ,λ) point 
 * acts as reference point (i.e. origin of ENU system).
 *
 * To perform the inverse operation (i.e. cartesian to topocentric) use the 
 * transpose of the R matrix, i.e. R^(T).
 *
 * @param[in] lat The geodetic or geocentric latitude in [rad]. Depending on 
 *            this vector, the corresponding u direction will be either 
 *            normal to the ellipsoid (at the given point), or normal to the 
 *            spheroid.
 * @param[in] lon Longitude of the reference point in [rad].
 * @return A 3x3 matrix; its first column is the e (East) unit vector, its 
 *         seconds column is the n (North) unit vector and the third column 
 *         is the u vector (Up). Hence:
 *         R = [e, n, u] where e, n, u are (3x1) unit vectors.
 */
Eigen::Matrix<double,3,3> geodetic2lvlh(double lat, double lon) noexcept;

}/* namespace dso */

#endif
