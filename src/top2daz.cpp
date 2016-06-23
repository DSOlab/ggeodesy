///
/// \file top2daz.cpp
///
/// \brief Compute azimouth, zenith and distance from a topocentric vector.
///

#include <cmath>
#include <stdexcept>
#include "geodesy.hpp"

void 
ngpt::top2daz(double north, double east, double up, double& distance, 
    double& azimouth, double& zenith)
{

  // spatial distance of vector
  distance  = std::sqrt(north*north + east*east + up*up);

  // check if zero distance or north are zero
  if ( (!distance) || (!north) ) {
      throw std::runtime_error("geodesy::top2daz -> Zero Division !!");
  }

  // azimouth
  double a { std::atan2(east, north) };

  // normalize to range [0-2pi)
  azimouth  = std::fmod(a, ngpt::D2PI);
  while (azimouth < .0) azimouth += ngpt::D2PI;

  // zenith angle [0-pi)
  zenith = std::acos(up / distance);

  // finished
  return;
}
