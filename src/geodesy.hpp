///
/// \file geodesy.hpp
///
/// \brief A list of frequently used geodetic functions.
///

#ifndef __GEODESY__
#define __GEODESY__

#include "geoconst.hpp"
#include <type_traits>

namespace ngpt
{

/// \brief Topocentric vector to azimouth, zenith and distance.
///
/// Compute the Distance, Azimouth and Zenith distance given a vector expressed
/// in a local, topocentric system (i.e. given the north, east and up 
/// components of the vector).
///
/// \param[in]  north    Vector north component, meters.
/// \param[in]  east     Vector east component, meters.
/// \param[in]  up       Vector up component, meters.
/// \param[out] distance The length of the vector, meters.
/// \param[out] azimouth The azimouth of the vector, radians [0,2*pi).
/// \param[out] zenith   The zenith distance, radians [0,pi)
/// \throw               std::runtime_error if zero division encountered.
///
/// \see "Physical Geodesy", pg. 210
///
void
top2daz(double north, double east, double up, double& distance,
    double& azimouth, double& zenith);

/// \brief Convert degrees to radians.
template<typename T> T deg2rad(T degrees) noexcept { return degrees * DEG2RAD; }

/// \brief Convert degrees to radians.
template<typename T> T rad2deg(T radians) noexcept { return radians * RAD2DEG; }

/// \brief Decimal to hexicondal degrees.
///
/// \tparam     T           Any floating point type for input and results.
/// \param[in]  decimal_deg The decimal degrees.
/// \param[out] deg         Integer degrees.
/// \param[out] min         Integer minutes.
/// \param[out] sec         The fractional seconds.
/// \throw                  Does not throw
///
/// \note In case a negative angle is given, then the (output) degrees are also
///       going to be negative.
///
template<typename T,
    typename = std::enable_if_t<
        std::is_floating_point<T>::value
        >
    >
    void decd2hexd(T decimal_deg, int& deg, int& min, T& sec) noexcept
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

/// \brief Hexicondal degrees to decimal degrees.
///
/// \tparam     T   Any floating point type for input and results.
/// \param[in] deg  Integer degrees.
/// \param[in] min  Integer minutes.
/// \param[in] sec  The fractional seconds.
/// \return         The angle in decimal degrees.
/// \throw          Does not throw
///
/// \note If the angle is negative, only the deg parameter should be negative;
///       if so the (decimal) degrees returned will also be negative. The other two
///       parameters are not checked for their sign, they should *ALWAYS* be
///       positive.
template<typename T,
    typename = std::enable_if_t<
        std::is_floating_point<T>::value
        >
    >
    T hexd2decd(int deg, int min, T sec) noexcept
{
    T angle { static_cast<T>(std::abs(deg)) +
        (static_cast<T>(min) + sec / 60.0e0) / 60.0e0 };

    return std::copysign(angle, (T)deg);
}

/// \brief Hexicondal degrees to radians.
///
/// \tparam     T   Any floating point type for input and results.
/// \param[in] deg  Integer degrees.
/// \param[in] min  Integer minutes.
/// \param[in] sec  The fractional seconds.
/// \return         The angle in radians.
/// \throw          Does not throw
///
/// \note If the angle is negative, only the deg parameter should be negative;
///       if so the radians returned will also be negative. The other two
///       parameters are not checked for their sign, they should *ALWAYS* be
///       positive.
template<typename T,
    typename = std::enable_if_t<
        std::is_floating_point<T>::value
        >
    >
    T hexd2rad(int deg, int min, T sec) noexcept
{
    return deg2rad( hexd2decd(deg, min, sec) );
}


/// \brief Radians to hexicondal degrees.
///
/// \tparam     T           Any floating point type for input and results.
/// \param[in]  radians     An angle in radians.
/// \param[out] deg         Integer degrees.
/// \param[out] min         Integer minutes.
/// \param[out] sec         The fractional seconds.
/// \throw                  Does not throw
///
/// \note In case a negative angle is given, then the (output) degrees are also
///       going to be negative.
///
template<typename T,
    typename = std::enable_if_t<
        std::is_floating_point<T>::value
        >
    >
    void rad2hexd(T radians, int& deg, int& min, T& sec) noexcept
{ return decd2hexd(rad2deg(radians), deg, min, sec); }

} // end namespace geodesy

#endif
