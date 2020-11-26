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
#include <boost/geometry.hpp>
#include <boost/geometry/formulas/karney_direct.hpp>

using namespace boost::geometry;
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
  ngpt::core::direct_karney(lat1, lon1, a1, s, A, F, B, lat2, lon2);

// For storing the resulting values.
formula::result_direct<double> result_drc;

// WGS-84 spheroid.
srs::spheroid<double> spheroid(A, B);

// Define the strategy.
typedef formula::karney_direct<double, true, true, false, false, 6u>
    karney_direct_type;
result_drc =
        karney_direct_type::apply(lon1, lat1, s, a1, spheroid);
//printf("\nLatitude  : %20.11f",rad2deg(result_drc.lat2));
//printf("\nLongtitude: %20.11f",rad2deg(result_drc.lon2));
//printf("\nAzimuth   : %20.11f",rad2deg(result_drc.reverse_azimuth));
printf("\nLatitude  : %20.11f",(result_drc.lat2));
printf("\nLongtitude: %20.11f",(result_drc.lon2));
printf("\nAzimuth   : %20.11f",(result_drc.reverse_azimuth));

  return 0;
}
