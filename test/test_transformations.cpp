#include "geodesy.hpp"
#include "test_help.hpp"
#include <cassert>
#include <iostream>
#include <random>

using dso::Ellipsoid;
using dso::ellipsoid;
using dso::VECTOR3;

constexpr auto wgs84 = Ellipsoid(ellipsoid::wgs84);

int main() {
  double lat1, lon1, hgt1, lat2, lon2, hgt2;
#ifdef CHECK_PRECISION
  double dlat, dlon, dhgt, dx, dy, dz;
  double max_dlat = std::numeric_limits<double>::min(),
         max_dlon = std::numeric_limits<double>::min(),
         max_dhgt = std::numeric_limits<double>::min(),
         max_dx = std::numeric_limits<double>::min(),
         max_dy = std::numeric_limits<double>::min(),
         max_dz = std::numeric_limits<double>::min();
#endif

    VECTOR3 ivec;
  // test ellipsoidal to cartesian and back
  for (int i = 0; i < 500; ++i) {
    lat1 = generate_random_double(-dso::DPI / 2e0, dso::DPI / 2e0);
    lon1 = generate_random_double(-dso::DPI, dso::DPI);
    hgt1 = generate_random_double(-10e0, 9e3);
    ivec(0) = lon1; ivec(1) = lat1; ivec(2) = hgt1;
    VECTOR3 v1 = dso::ell2car<ellipsoid::wgs84>(ivec);
    VECTOR3 v2 = dso::ell2car(ivec, wgs84);
#ifdef CHECK_PRECISION
    if ((dx = std::abs(v1(0) - v2(0))) > max_dx)
      max_dx = dx;
    if ((dy = std::abs(v1(1) - v2(1))) > max_dy)
      max_dy = dy;
    if ((dz = std::abs(v1(2) - v2(2))) > max_dz)
      max_dz = dz;
#else
    assert(approxEqual(v1(0), v2(0)) && approxEqual(v1(1), v2(1)) && approxEqual(v1(2), v2(2)));
#endif
    
    v2=dso::car2ell<ellipsoid::wgs84>(v1);
    lon2 = v2(0); lat2=v2(1); hgt2=v2(2);
#ifdef CHECK_PRECISION
    if ((dlat = std::abs(m_rad2meters<ellipsoid::wgs84>(lat1 - lat2, lat1))) >
        max_dlat)
      max_dlat = dlat;
    if ((dlon = std::abs(p_rad2meters<ellipsoid::wgs84>(lon1 - lon2, lat1))) >
        max_dlon)
      max_dlon = dlon;
    if ((dhgt = std::abs(hgt1 - hgt2)) > max_dhgt)
      max_dhgt = dhgt;
#else
    assert(approxEqual(lat1, lat2) && approxEqual(lon1, lon2) &&
           approxEqual(hgt1, hgt2));
#endif
    // dso::car2ell(x1,y1,z1,wgs84,lat1,lon1,hgt1);
    // assert(approxEqual(lat1,lat2) && approxEqual(lon1,lon2) &&
    // approxEqual(hgt1,hgt2));
  }

#ifdef CHECK_PRECISION
  printf("\nMax values for error:");
  printf("\n\tMax dx   = %20.15fm", max_dx);
  printf("\n\tMax dy   = %20.15fm", max_dy);
  printf("\n\tMax dz   = %20.15fm", max_dz);
  printf("\n\tMax dlat = %20.15fm", max_dlat);
  printf("\n\tMax dlon = %20.15fm", max_dlon);
  printf("\n\tMax dhgt = %20.15fm", max_dhgt);
#endif

  return 0;
}
