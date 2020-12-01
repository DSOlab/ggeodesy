#include "geodesy.hpp"
#include "test_help.hpp"
#include "vincenty.hpp"
#include "units.hpp"
#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <iostream>
#include <random>

using namespace ngpt;

int main() {

constexpr double A = ellipsoid_traits<ellipsoid::wgs84>::a,
                 F = ellipsoid_traits<ellipsoid::wgs84>::f,
                 B = semi_minor<ellipsoid::wgs84>();

  const double lat1=40e0;
  const double lon1=0e0;
  const double a1=30e0;
  const double s=10000000e0;
  double lat2, lon2;
  double a2 = ngpt::core::direct_karney(lat1, lon1, a1, s, A, F, B, lat2, lon2);

//printf("\nLatitude  : %20.11f",rad2deg(result_drc.lat2));
//printf("\nLongtitude: %20.11f",rad2deg(result_drc.lon2));
//printf("\nAzimuth   : %20.11f",rad2deg(result_drc.reverse_azimuth));
printf("\nLatitude  : %20.11f",(lat2));
printf("\nLongtitude: %20.11f",(lon2));
printf("\nAzimuth   : %20.11f",(a2));

  return 0;
}
