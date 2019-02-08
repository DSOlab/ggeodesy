///
/// @file geodesy.hpp
///
/// @brief A list of frequently used geodetic functions.
///
/// @see http://www.movable-type.co.uk/scripts/latlong.html

#ifndef __NGPT_GEODESY_HPP__
#define __NGPT_GEODESY_HPP__

#include <type_traits>
#include <cmath>
#include <stdexcept>
#include "geoconst.hpp"

namespace ngpt
{

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
void
top2daz(double north, double east, double up, double& distance,
    double& azimouth, double& zenith);

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
void
top2car(double north, double east, double up, double lat, double lon,
    double& dx, double& dy, double& dz) noexcept;

/// @brief Convert degrees to radians.
/// @tparam    T       Any floating type
/// @param[in] degrees Angle in decimal degrees
/// @return            The (input) angle in radians
template<typename T>
    T
    deg2rad(T degrees) noexcept
    { return degrees * DEG2RAD; }

/// @brief Convert radians to degrees.
/// @tparam    T       Any floating type
/// @param[in] radians Angle in radians
/// @return            The (input) angle in decimal degrees
template<typename T>
    T
    rad2deg(T radians) noexcept
    { return radians * RAD2DEG; }

/// @brief Normalize angle.
///
/// Normalize an angle in the interval [lower, upper). 
///
/// @tparam    T      Any floating point type for input and results.
/// @param[in] angle  The angle to normalize (note that the unit should be the
///                   same as in lower and upper parameters).
/// @param[in] lower  lower bound (inclusive). Default value is 0
/// @param[in] upper  upper bound (exclusive). Default is 2*π
///
/// @note  It is not always needed to use this function to normalize an angle.
///        If i.e. the angle is a result of a function that returns values in
///        the range (-π , +π] radians, and we need to normalize in the range
///        [0, 2π), then we can use: angle = std::fmod(angle+2π, 2π).
template<typename T,
    typename = std::enable_if_t<
        std::is_floating_point<T>::value
        >
    >
    T
    normalize_angle(T angle, T lower=0e0, T upper=D2PI)
{
    if (lower >= upper)
        throw std::invalid_argument("[ERROR] normalize_angle(): Invalid lower/upper bounds");

    double res {angle};
    if (angle>upper || angle==lower)
        angle = lower +
            std::fmod(std::abs(angle+upper), std::abs(lower)+std::abs(upper));
    if (angle<lower || angle==upper)
        angle = upper -
            std::fmod(std::abs(angle-lower), std::abs(lower)+std::abs(upper));

    res = (res==upper)?(lower):(angle);
    return res;
}

/// @brief Decimal to hexicondal degrees.
///
/// @tparam     T           Any floating point type for input and results.
/// @param[in]  decimal_deg The decimal degrees.
/// @param[out] deg         Integer degrees.
/// @param[out] min         Integer minutes.
/// @param[out] sec         The fractional seconds.
/// @throw                  Does not throw
///
/// @note In case a negative angle is given, then the (output) degrees are also
///       going to be negative.
///
template<typename T,
    typename = std::enable_if_t<
        std::is_floating_point<T>::value
        >
    >
    void
    decd2hexd(T decimal_deg, int& deg, int& min, T& sec) noexcept
{
    T decdeg { std::abs(decimal_deg) };

    deg = static_cast<int>(decdeg);
    min = static_cast<int>( (decdeg - static_cast<T>(deg)) *
        static_cast<T>(60.0e0) );
    sec = decdeg - (static_cast<T>(deg) + static_cast<T>(min)/60.0e0);
    sec *= 3600.0e0;

    if (decimal_deg < T{0}) deg *= -1;
    return;
}

/// @brief Hexicondal degrees to decimal degrees.
///
/// @tparam     T   Any floating point type for input and results.
/// @param[in] deg  Integer degrees.
/// @param[in] min  Integer minutes.
/// @param[in] sec  The fractional seconds.
/// @return         The angle in decimal degrees.
/// @throw          Does not throw
///
/// @note If the angle is negative, only the deg parameter should be negative;
///       if so the (decimal) degrees returned will also be negative. The other two
///       parameters are not checked for their sign, they should *ALWAYS* be
///       positive.
template<typename T,
    typename = std::enable_if_t<
        std::is_floating_point<T>::value
        >
    >
    T
    hexd2decd(int deg, int min, T sec) noexcept
{
    T angle { static_cast<T>(std::abs(deg)) +
        (static_cast<T>(min) + sec / 60.0e0) / 60.0e0 };

    return std::copysign(angle, (T)deg);
}

/// @brief Hexicondal degrees to radians.
///
/// @tparam     T   Any floating point type for input and results.
/// @param[in] deg  Integer degrees.
/// @param[in] min  Integer minutes.
/// @param[in] sec  The fractional seconds.
/// @return         The angle in radians.
/// @throw          Does not throw
///
/// @note If the angle is negative, only the deg parameter should be negative;
///       if so the radians returned will also be negative. The other two
///       parameters are not checked for their sign, they should *ALWAYS* be
///       positive.
template<typename T,
    typename = std::enable_if_t<
        std::is_floating_point<T>::value
        >
    >
    T
    hexd2rad(int deg, int min, T sec) noexcept
{return deg2rad( hexd2decd(deg, min, sec) );}


/// @brief Radians to hexicondal degrees.
///
/// @tparam     T           Any floating point type for input and results.
/// @param[in]  radians     An angle in radians.
/// @param[out] deg         Integer degrees.
/// @param[out] min         Integer minutes.
/// @param[out] sec         The fractional seconds.
/// @throw                  Does not throw
///
/// @note In case a negative angle is given, then the (output) degrees are also
///       going to be negative.
///
template<typename T,
    typename = std::enable_if_t<
        std::is_floating_point<T>::value
        >
    >
    void
    rad2hexd(T radians, int& deg, int& min, T& sec) noexcept
{return decd2hexd(rad2deg(radians), deg, min, sec);}

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
template<typename T,
    typename = std::enable_if_t<
        std::is_floating_point<T>::value
        >
    >
    T
    bearing(T lat1, T lon1, T lat2, T lon2) noexcept
{
    using std::sin;
    using std::cos;

    double DeltaLambda = lon2 - lon1;
    double cosLat2     = cos(lat2);
    double nom         = sin(DeltaLambda) * cosLat2;
    double denom       = cos(lat1)*sin(lat2) -
                         sin(lat1)*cosLat2*cos(DeltaLambda);
    double angle       = std::atan2(nom, denom);
    // normalize to [0-2pi)
    return std::fmod(angle+D2PI, D2PI);
}

} // end namespace geodesy

#endif
