/// @file top2car.cpp
#include "geodesy.hpp"
#include <cmath>

using dso::VECTOR3;

void dso::top2car(double east, double north, double up, double lon, double lat,
                  double &dx, double &dy, double &dz) noexcept {
  const double slon = std::sin(lon);
  const double clon = std::cos(lon);
  const double slat = std::sin(lat);
  const double clat = std::cos(lat);

  dx = -slon * east - clon * slat * north + clon * clat * up;
  dy = clon * east - slon * slat * north + slon * clat * up;
  dz = clat * north + slat * up;

  return;
}
