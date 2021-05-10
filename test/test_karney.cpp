#include "geodesy.hpp"
#include "test_help.hpp"
#include "units.hpp"
#include "vincenty.hpp"
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

  double lat1 = ngpt::deg2rad(40e0);
  double lon1 = 0e0;
  double a1 = ngpt::deg2rad(30e0);
  double s = 10000000e0;
  double lat2, lon2;
  double a2 = ngpt::core::direct_karney(lat1, lon1, a1, s, A, F, B, lat2, lon2);

  // printf("\nLatitude  : %20.11f",rad2deg(result_drc.lat2));
  // printf("\nLongtitude: %20.11f",rad2deg(result_drc.lon2));
  // printf("\nAzimuth   : %20.11f",rad2deg(result_drc.reverse_azimuth));
  printf("\nLatitude  : %20.11f", ngpt::rad2deg(lat2));
  printf("\nLongtitude: %20.11f", ngpt::rad2deg(lon2));
  printf("\nAzimuth   : %20.11f", ngpt::rad2deg(a2));

  lat1 = ngpt::deg2rad(-30e0);
  lat2 = ngpt::deg2rad(29.9e0);
  lon2 = 0e0;
  lon1 = ngpt::deg2rad(179.8e0);
  double a12, a21;
  ngpt::core::inverse_karney(lat1, lon1, lat2, lon2,  A, F, B, a12,  a21);

  return 0;
}
