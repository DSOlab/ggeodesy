#include "geodesy.hpp"
#include "units.hpp"

using dso::VECTOR3;

void dso::top2dae(const VECTOR3 &enu, double &distance, double &azimouth,
                  double &elevation) {
  const double e = enu(0);
  const double n = enu(1);
  const double u = enu(2);

  const double rho2 = e * e + n * n;
  const double rho = std::sqrt(rho2);

  // azimouth in [0,2π]
  azimouth = norm_angle<double, dso::AngleUnit::Radians>(std::atan2(e, n));

  // elevation angle [0-π]
  elevation = std::atan(u / rho);

  // distance
  distance = rho;
  return;
}

void dso::top2dae(const VECTOR3 &enu, double &distance, double &azimouth,
                  double &elevation, VECTOR3 &dAdr, VECTOR3 &dEdr) {
  const double e = enu(0);
  const double n = enu(1);
  const double u = enu(2);

  const double rho2 = e * e + n * n;
  const double rho = std::sqrt(rho2);

  // azimouth in [0,2π]
  azimouth = norm_angle<double, dso::AngleUnit::Radians>(std::atan2(e, n));

  // elevation angle [0-π]
  elevation = std::atan(u / rho);

  // partials
  dAdr(0) = n / rho2;
  dAdr(1) = -e / rho2;
  dAdr(2) = 0e0;
  dEdr(0) = -e * u / rho2 / rho;
  dEdr(1) = -n * u / rho2 / rho;
  dEdr(2) = 1e0 / rho;

  // distance
  distance = rho;
  return;
}
