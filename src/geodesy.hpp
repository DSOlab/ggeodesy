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

/// \brief Decimal to hexicondal degrees.
template<typename T,
        std::enable_if<
            std::is_floating_point<T>::value,
            true
            >
        >
    decd2hexd(

} // end namespace geodesy

#endif
