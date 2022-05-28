/// @file top2daz.cpp
///
/// @brief Compute azimouth, zenith and distance from a topocentric vector.
///
/// @see
/// http://www.navipedia.net/index.php/Transformations_between_ECEF_and_ENU_coordinates

#include "geodesy.hpp"
#include "units.hpp"
#include <cmath>
#include <stdexcept>

void dso::top2daz(double east, double north, double up, double &distance,
                   double &azimouth, double &zenith) {

  // spatial distance of vector
  distance = std::sqrt(north * north + east * east + up * up);

  // check if zero distance or north are zero
  if ((!distance) || (!north)) {
    throw std::runtime_error("[ERROR] geodesy::top2daz -> Zero Division !!");
  }

  // azimouth in [0,2π]
  azimouth = norm_angle<double, AngleUnit::Radians>(std::atan2(east, north));

  // zenith angle [0-π]
  zenith = std::acos(up / distance);

  // finished
  return;
}

void dso::top2dae(const double *enu, double &distance, double &azimouth,
                  double &elevation) {
  const double e = enu[0];
  const double n = enu[1];
  const double u = enu[2];

  const double rho2 = e*e + n*n;
  const double rho = std::sqrt(rho2);

  // azimouth in [0,2π]
  azimouth = norm_angle<double, AngleUnit::Radians>(std::atan2(e, n));
  
  // elevation angle [0-π]
  elevation = std::atan(u/rho);

  // distance
  distance = rho;
  return;
}

void dso::top2dae(const double *enu, double &distance, double &azimouth,
                  double &elevation, double *dAdr, double *dEdr) {
  const double e = enu[0];
  const double n = enu[1];
  const double u = enu[2];

  const double rho2 = e*e + n*n;
  const double rho = std::sqrt(rho2);

  // azimouth in [0,2π]
  azimouth = norm_angle<double, AngleUnit::Radians>(std::atan2(e, n));
  
  // elevation angle [0-π]
  elevation = std::atan(u/rho);

  // partials
  dAdr[0]=n/rho2; dAdr[1] = -e/rho2; dAdr[2]=0e0;
  dEdr[0]=-e*u/rho2/rho; dEdr[1] = -n*u/rho2/rho; dEdr[2]=1e0/rho;

  // distance
  distance = rho;
  return;
}
