#include "geoconst.hpp"
#include "geodesy.hpp"
#include "units.hpp"
#include "vincenty.hpp"
#include <cmath>
#include <stdexcept>

using dso::D2PI;
using dso::DPI;

/// @brief arrange lat1 and lat2 such that:
///        phi1 <= 0 and phi1 <= phi2 <= -phi1
///        This is the "canonical configuration" as described in Karney 2013,
///        eq. 44
void arrange_lats(double lat1, double lat2, double &phi1, double& phi2) noexcept {
  if (lat1>=0e0 && lat2>=0e0) {
    phi1 = -std::max(lat1, lat2);
    phi2 = std::min(lat1, lat2);
  } else if (lat1<=0e0 && lat2>0e0) {
    if (std::abs(lat1)>std::abs(lat2)) {
      phi1 = lat1;
      phi2 = lat2;
    } else {
      phi1 = -lat2;
      phi2 = lat1;
    }
  } else if (lat1>=0e0 && lat2<0e0) {
    if (std::abs(lat1)>std::abs(lat2)) {
      phi1 = -lat1;
      phi2 = lat2;
    } else {
      phi1 = lat2;
      phi2 = lat1;
    }
  } else {
    phi1 = -std::max(-lat1, -lat2);
    phi2 = -std::min(-lat1, -lat2);
  }
#ifdef DEBUG
  assert(phi1<=0e0);
  assert(phi1<=phi2 && phi2<=-phi1);
#endif
  return;
}

void arrange_lats_boost(double lat1, double lat2, double &phi1, double& phi2) noexcept {
  int swap_point = std::abs(lat1) < std::abs(lat2) ? -1 : 1;
  if (swap_point < 0) {
    std::swap(lat1,lat2);
  }
  int lat_sign = lat1 < 0 ? 1 : -1;
  lat1 *= lat_sign;
  lat2 *= lat_sign;
  phi1 = lat1;
  phi2 = lat2;
  return;
}

double dso::core::inverse_karney(double lat1, double lon1, double lat2,
                                    double lon2, double semi_major,
                                    double flattening, double semi_minor,
                                    double &a12, double &a21
                                    ) {
  const double a = semi_major;
  const double f = flattening;
  const double b = semi_minor;
  assert(b!=0e0); // delete this shit

  double phi1, phi2;
  // arrange latitude in canonical configuration, eq. 44
  arrange_lats(lat1, lat2, phi1, phi2);
#ifdef DEBUG
  double p1,p2;
  arrange_lats_boost(lat1, lat2, p1, p2);
  assert(phi1==p1 && phi2==p2);
#endif
  const double l12 = normalize_angle(lon2-lon1, 0e0, DPI);
#ifdef DEBUG
  assert(l12>=0e0 && l12 <= DPI);
#endif

  // solve the astroid problem
  const double beta1 = dso::core::reduced_latitude(lat1, f);
  const double beta2 = dso::core::reduced_latitude(lat2, f);
  const double cosbeta1 = std::cos(beta1);
  const double delta = f*a*DPI*cosbeta1*cosbeta1;
  const double x = (l12 - DPI)*(a*cosbeta1)/delta;
  printf("\nx=%+10.6f", x);
  const double y = (beta2 + beta1)*a/delta;
  printf("\ny=%+10.6f", y);

  a12 = a21 = 0e0;
  return 0e0;
}
