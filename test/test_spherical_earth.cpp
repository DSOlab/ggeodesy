#include "geodesy.hpp"
#include "great_circle.hpp"
#include "test_help.hpp"
#include <algorithm>
#include <array>
#include <cassert>
#include <iostream>
#include <random>

int main() {

  double lat1, lat2, lon1, lon2, d1, d2, d3, d4, d5, d6, meand, dd1, dd2, dd3,
      dd4, dd5, dd6;
  const double R = ngpt::mean_earth_radius<ngpt::ellipsoid::wgs84>();
  std::vector<std::array<double, 7>> diffs;

  for (int i = 0; i < 50; ++i) {
    // random point 1; make point 2 too close
    lat1 = generate_random_double(-ngpt::DPI / 2e0, ngpt::DPI / 2e0);
    lat2 = lat1 + generate_random_double(-1e-3, 1e-3);
    lon1 = generate_random_double(-ngpt::DPI, ngpt::DPI);
    lon2 = lon1 + generate_random_double(-1e-3, 1e-3);
    // compute great circle distance
    d1 = ngpt::core::great_circle_distance_cosines(lat1, lon1, lat2, lon2, R);
    d2 = ngpt::core::great_circle_distance_haversine(lat1, lon1, lat2, lon2, R);
    d3 = ngpt::core::great_circle_distance_vincenty(lat1, lon1, lat2, lon2, R);
    d4 = ngpt::core::great_circle_distance_chordl(lat1, lon1, lat2, lon2, R);
    d5 =
        ngpt::core::great_circle_distance_pythagoras(lat1, lon1, lat2, lon2, R);
    d6 = ngpt::core::great_circle_distance_polar(lat1, lon1, lat2, lon2, R);
    // printf("\n%15.5f %15.5f %15.5f %15.5f %15.5f %15.5f", d1, d2, d3, d4, d5,
    // d6);
    meand = (d1 + d2 + d3 + d4) / 4e0;
    dd1 = std::abs(d1 - meand);
    dd2 = std::abs(d2 - meand);
    dd3 = std::abs(d3 - meand);
    dd4 = std::abs(d4 - meand);
    dd5 = std::abs(d5 - meand);
    dd6 = std::abs(d6 - meand);
    // printf(" %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f", dd1, dd2, dd3, dd4,
    // dd5, dd6);
    diffs.emplace_back(
        std::array<double, 7>{meand, dd1, dd2, dd3, dd4, dd5, dd6});
  }

  for (int i = 0; i < 50; ++i) {
    // random point 1; make point 2 too close
    lat1 = generate_random_double(-ngpt::DPI / 2e0, ngpt::DPI / 2e0);
    lat2 = lat1 + generate_random_double(-1e-2, 1e-2);
    lon1 = generate_random_double(-ngpt::DPI, ngpt::DPI);
    lon2 = lon1 + generate_random_double(-1e-2, 1e-2);
    // compute great circle distance
    d1 = ngpt::core::great_circle_distance_cosines(lat1, lon1, lat2, lon2, R);
    d2 = ngpt::core::great_circle_distance_haversine(lat1, lon1, lat2, lon2, R);
    d3 = ngpt::core::great_circle_distance_vincenty(lat1, lon1, lat2, lon2, R);
    d4 = ngpt::core::great_circle_distance_chordl(lat1, lon1, lat2, lon2, R);
    d5 =
        ngpt::core::great_circle_distance_pythagoras(lat1, lon1, lat2, lon2, R);
    d6 = ngpt::core::great_circle_distance_polar(lat1, lon1, lat2, lon2, R);
    // printf("\n%15.5f %15.5f %15.5f %15.5f %15.5f %15.5f", d1, d2, d3, d4, d5,
    // d6);
    meand = (d1 + d2 + d3 + d4) / 4e0;
    dd1 = std::abs(d1 - meand);
    dd2 = std::abs(d2 - meand);
    dd3 = std::abs(d3 - meand);
    dd4 = std::abs(d4 - meand);
    dd5 = std::abs(d5 - meand);
    dd6 = std::abs(d6 - meand);
    // printf(" %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f", dd1, dd2, dd3, dd4,
    // dd5, dd6);
    diffs.emplace_back(
        std::array<double, 7>{meand, dd1, dd2, dd3, dd4, dd5, dd6});
  }

  for (int i = 0; i < 50; ++i) {
    // random point 1; make point 2 too close
    lat1 = generate_random_double(-ngpt::DPI / 2e0, ngpt::DPI / 2e0);
    lat2 = lat1 + generate_random_double(-1e-1, 1e-1);
    lon1 = generate_random_double(-ngpt::DPI, ngpt::DPI);
    lon2 = lon1 + generate_random_double(-1e-1, 1e-1);
    // compute great circle distance
    d1 = ngpt::core::great_circle_distance_cosines(lat1, lon1, lat2, lon2, R);
    d2 = ngpt::core::great_circle_distance_haversine(lat1, lon1, lat2, lon2, R);
    d3 = ngpt::core::great_circle_distance_vincenty(lat1, lon1, lat2, lon2, R);
    d4 = ngpt::core::great_circle_distance_chordl(lat1, lon1, lat2, lon2, R);
    d5 =
        ngpt::core::great_circle_distance_pythagoras(lat1, lon1, lat2, lon2, R);
    d6 = ngpt::core::great_circle_distance_polar(lat1, lon1, lat2, lon2, R);
    // printf("\n%15.5f %15.5f %15.5f %15.5f %15.5f %15.5f", d1, d2, d3, d4, d5,
    // d6);
    meand = (d1 + d2 + d3 + d4) / 4e0;
    dd1 = std::abs(d1 - meand);
    dd2 = std::abs(d2 - meand);
    dd3 = std::abs(d3 - meand);
    dd4 = std::abs(d4 - meand);
    dd5 = std::abs(d5 - meand);
    dd6 = std::abs(d6 - meand);
    // printf(" %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f", dd1, dd2, dd3, dd4,
    // dd5, dd6);
    diffs.emplace_back(
        std::array<double, 7>{meand, dd1, dd2, dd3, dd4, dd5, dd6});
  }

  for (int i = 0; i < 1500; ++i) {
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
    d5 =
        ngpt::core::great_circle_distance_pythagoras(lat1, lon1, lat2, lon2, R);
    d6 = ngpt::core::great_circle_distance_polar(lat1, lon1, lat2, lon2, R);
    // printf("\n%15.5f %15.5f %15.5f %15.5f %15.5f %15.5f", d1, d2, d3, d4, d5,
    // d6);
    meand = (d1 + d2 + d3 + d4) / 4e0;
    dd1 = std::abs(d1 - meand);
    dd2 = std::abs(d2 - meand);
    dd3 = std::abs(d3 - meand);
    dd4 = std::abs(d4 - meand);
    dd5 = std::abs(d5 - meand);
    dd6 = std::abs(d6 - meand);
    // printf(" %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f", dd1, dd2, dd3, dd4,
    // dd5, dd6);
    diffs.emplace_back(
        std::array<double, 7>{meand, dd1, dd2, dd3, dd4, dd5, dd6});
  }

  printf("\n%12s %12s %12s %12s %12s %12s %12s", "Av. Distance", "Cosines",
         "Haversine", "Vincenty", "ChordL", "Pythagoras", "Polar");
  std::sort(diffs.begin(), diffs.end(),
            [](const std::array<double, 7> &a1,
               const std::array<double, 7> &a2) { return a1[0] < a2[0]; });
  for (const auto &a : diffs)
    printf("\n%12.3f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f", a[0] / 1e3,
           a[1], a[2], a[3], a[4], a[5], a[6]);

  printf("\n");
  return 0;
}
