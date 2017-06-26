///
/// \file top2car.cpp
///

#include <cmath>
#include "geodesy.hpp"

void 
ngpt::top2car(double north, double east, double up, double lat, double lon,
    double& dx, double& dy, double& dz)
noexcept
{
    const double slon = std::sin(lon);
    const double clon = std::cos(lon);
    const double slat = std::sin(lat);
    const double clat = std::cos(lat);

    dx = -slon*east -clon*slat*north + clon*clat*up;
    dy =  clon*east -slon*slat*north + slon*clat*up;
    dz =             clat*north      + slat*up;

    return;
}
