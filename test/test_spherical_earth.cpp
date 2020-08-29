#include "geodesy.hpp"
#include "great_circle.hpp"
#include "test_help.hpp"
#include <cassert>
#include <iostream>
#include <random>

int main() {

  double lat1, lat2, lon1, lon2, d1, d2, d3, d4, meand, dd1, dd2, dd3, dd4;
  const double R = ngpt::mean_earth_radius<ngpt::ellipsoid::wgs84>();

  for (int i = 0; i < 150; ++i) {
    // random points 1 and 2
    lat1 = generate_random_double(-ngpt::DPI / 2e0, ngpt::DPI / 2e0);
    lat2 = generate_random_double(-ngpt::DPI / 2e0, ngpt::DPI / 2e0);
    lon1 = generate_random_double(-ngpt::DPI, ngpt::DPI);
    lon2 = generate_random_double(-ngpt::DPI, ngpt::DPI);
    // compute great circle distance
    d1 = ngpt::core::great_circle_distance_cosines(lat1, lon1, lat2, lon2, R);
    d2 = ngpt::core::great_circle_distance_haversine(lat1, lon1, lat2, lon2, R);
    d3 = ngpt::core::great_circle_distance_vincenty(lat1, lon1, lat2, lon2, R);
    d4 = ngpt::core::great_circle_distance_chordl(lat1, lon1, lat2, lon2, R);
    printf("\n%15.5f %15.5f %15.5f %15.5f", d1, d2, d3, d4);
    meand = (d1 + d2 + d3 + d4) / 4e0;
    dd1 = std::abs(d1 - meand);
    dd2 = std::abs(d2 - meand);
    dd3 = std::abs(d3 - meand);
    dd4 = std::abs(d4 - meand);
    printf(" %12.8f %12.8f %12.8f %12.8f", dd1, dd2, dd3, dd4);
  }

  printf("\n");
  return 0;
}
