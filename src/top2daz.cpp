/// @file top2daz.cpp
/// @brief Compute azimouth, zenith and distance from a topocentric vector.
/// @see
/// http://www.navipedia.net/index.php/Transformations_between_ECEF_and_ENU_coordinates
#include "geodesy.hpp"
#include "units.hpp"

double dso::top2daz(double east, double north, double up,
                  double &azimouth, double &zenith) {

  // spatial distance of vector
  const double distance = std::sqrt(north * north + east * east + up * up);

  // check if zero distance or north are zero
  if ((!distance) || (!north)) {
    throw std::runtime_error("[ERROR] geodesy::top2daz -> Zero Division !!");
  }

  // azimouth in [0,2π]
  azimouth = norm_angle<double, AngleUnit::Radians>(std::atan2(east, north));

  // zenith angle [0-π]
  zenith = std::acos(up / distance);

  // finished
  return distance;
}
