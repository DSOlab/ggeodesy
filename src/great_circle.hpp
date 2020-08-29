///
/// @file great_circle.hpp
///
/// @brief Algorithms to compute great-circle distance or orthodromic distance
///        between two points on the surface of a sphere, measured along the
///        surface of the sphere.
///        The Earth is nearly spherical, so great-circle distance formulas
///        give the distance between points on the surface of the Earth correct
///        to within about 0.5%.
///        Let λ1, ϕ1 and λ2, ϕ2 be the geographical longitude and latitude in
///        radians of two points 1 and 2, and Δλ, Δϕ be their absolute
///        differences; the key thing then for arc legth computation is to
///        find Δσ, the central angle between them. Once we know Δσ, the actual
///        arc length d on a sphere of radius r can be trivially computed as
///        d = r * Δσ
///
/// @see   https://en.wikipedia.org/wiki/Great-circle_distance
/// @todo

#ifndef __NGPT_SPHERICAL_EARTH_HPP__
#define __NGPT_SPHERICAL_EARTH_HPP__

#include <cmath>
#include <stdexcept>

namespace ngpt {

namespace core {
/// @brief Compute the haversine function.
///
/// This function is used within the haversine formula to compute great-circle
/// distances between two points on a sphere from their longitudes and
/// latitudes.
/// It is called the "haversine function" (half a versine) and given an angle
/// θ it is computed as:
/// \f$ hav(\theta ) =  \sin^2 \frac{\theta}{2} = \frac{1- \cos \theta}{2}\f$
///
/// @param[in] angle  The angle to compute the haversine function for (radians)
/// @return The haversine function
///
/// @see https://en.wikipedia.org/wiki/Haversine_formula
double haversine_angle(double angle) noexcept;

/// @brief Compute great circle distance between two points on a sphere of
///        radius R
/// @param[in] lat1 Latitude of starting point (radians)
/// @param[in] lon1 Longtitude of starting point (radians)
/// @param[in] lat2 Latitude of ending point (radians)
/// @param[in] lon2 Longtitude of ending point (radians)
/// @param[in] R    Radius of sphere (meters)
/// @return         Great circle distance between points 1 and 2 (meters)
///
/// This implementation uses the spherical law of cosines to compute the arc
/// length, using the central angle between them.
/// On computer systems with low floating-point precision, the spherical law
/// of cosines formula can have large rounding errors if the distance is small
/// (if the two points are a kilometer apart on the surface of the Earth, the
/// cosine of the central angle is near 0.99999999). For modern 64-bit
/// floating-point numbers, the spherical law of cosines formula, does not
/// have serious rounding errors for distances larger than a few meters on the
/// surface of the Earth. The haversine formula is numerically
/// better-conditioned for small distances. see
/// https://en.wikipedia.org/wiki/Great-circle_distance
double great_circle_distance_cosines(double lat1, double lon1, double lat2,
                                     double lon2, double R) noexcept;

/// @brief Compute great circle distance between two points on a sphere of
///        radius R
/// @param[in] lat1  Latitude of point 1 (radians)
/// @param[in] lon1  Longtitude  point 1 (radians)
/// @param[in] lat2  Latitude of point 2 (radians)
/// @param[in] lon2  Longtitude  point 2 (radians)
/// @param[in] R     Radius of sphere (meters)
/// @return          Great circle distance between points 1 and 2 (meters)
///
/// This implementation uses the haversine function to compute the arc length
/// between two points on a sphere. Although this formula is accurate for most
/// distances on a sphere, it too suffers from rounding errors for the special
/// (and somewhat unusual) case of antipodal points (on opposite ends of the
/// sphere).
///
/// @warning h only approaches 1 for antipodal points (on opposite sides of the
///          sphere). This will cause the formula to fail!
///
/// @note The haversine formula and law of cosines can't be guaranteed correct
///       to better than 0.5%
///
/// @see https://en.wikipedia.org/wiki/Haversine_formula
/// https://en.wikipedia.org/wiki/Great-circle_distance
double great_circle_distance_haversine(double lat1, double lon1, double lat2,
                                       double lon2, double R) noexcept;

/// @brief Compute great circle distance between two points on a sphere of
///        radius R
/// @param[in] lat1  Latitude of point 1 (radians)
/// @param[in] lon1  Longtitude  point 1 (radians)
/// @param[in] lat2  Latitude of point 2 (radians)
/// @param[in] lon2  Longtitude  point 2 (radians)
/// @param[in] R     Radius of sphere (meters)
/// @return          Great circle distance between points 1 and 2 (meters)
///
/// This implementation uses a special case of the Vincenty formula for an
/// ellipsoid with equal major and minor axes to compute the central angle
/// (between the two points). formula This formula is accurate for all distances
/// including antipodal points.
/// https://en.wikipedia.org/wiki/Great-circle_distance
double great_circle_distance_vincenty(double lat1, double lon1, double lat2,
                                      double lon2, double R) noexcept;

/// @brief Compute great circle distance between two points on a sphere of
///        radius R
/// @param[in] lat1  Latitude of point 1 (radians)
/// @param[in] lon1  Longtitude  point 1 (radians)
/// @param[in] lat2  Latitude of point 2 (radians)
/// @param[in] lon2  Longtitude  point 2 (radians)
/// @param[in] R     Radius of sphere (meters)
/// @return          Great circle distance between points 1 and 2 (meters)
///
/// This implementation uses the great circle chord length calculated for the
/// corresponding unit sphere, by means of Cartesian subtraction
/// https://en.wikipedia.org/wiki/Great-circle_distance
double great_circle_distance_chordl(double lat1, double lon1, double lat2,
                                    double lon2, double R) noexcept;

} // namespace core

} // namespace ngpt

#endif
