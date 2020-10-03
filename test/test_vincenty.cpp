#include "geodesy.hpp"
#include "test_help.hpp"
#include "vincenty.hpp"
#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <iostream>
#include <random>

using namespace ngpt;

double rad2seconds(double rad) noexcept {
  double ddeg = rad2deg(rad);
  return ddeg * 3600e0;
}

int main() {

  double lat1, lon1, s, az, lat21, lon21, lat22, lon22, baz1, baz2;
  const double R = ngpt::mean_earth_radius<ngpt::ellipsoid::wgs84>();
  const double a = ellipsoid_traits<ellipsoid::wgs84>::a,
               f = ellipsoid_traits<ellipsoid::wgs84>::f,
               b = semi_minor<ellipsoid::wgs84>();

  for (int i = 0; i < 1500; ++i) {
    // random points 1 and 2
    lat1 = generate_random_double(-ngpt::DPI / 2e0, ngpt::DPI / 2e0);
    // lat2 = generate_random_double(-ngpt::DPI / 2e0, ngpt::DPI / 2e0);
    lon1 = generate_random_double(-ngpt::DPI, ngpt::DPI);
    // lon2 = generate_random_double(-ngpt::DPI, ngpt::DPI);
    s = generate_random_double(0e0, ngpt::D2PI * R);
    az = generate_random_double(0e0, ngpt::D2PI);

    baz1 =
        ngpt::core::direct_vincenty(lat1, lon1, az, s, a, f, b, lat21, lon21);
    baz2 =
        ngpt::core::direct_vincenty2(lat1, lon1, az, s, a, f, b, lat22, lon22);
    printf("\n%+20.15f %+20.15f %+20.15f %+15.5f %+15.5f",
           rad2seconds(baz1 - baz2), rad2seconds(lat21 - lat22),
           rad2seconds(lon21 - lon22),
           infinitesimal_meridian_arc<ellipsoid::wgs84, double>(lat21,
                                                                lat21 - lat22),
           parallel_arc_length<ellipsoid::wgs84, double>(lat21, lon21 - lon22));
  }

  auto t1 = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < 15000; ++i) {
    // random points 1 and 2
    lat1 = generate_random_double(-ngpt::DPI / 2e0, ngpt::DPI / 2e0);
    // lat2 = generate_random_double(-ngpt::DPI / 2e0, ngpt::DPI / 2e0);
    lon1 = generate_random_double(-ngpt::DPI, ngpt::DPI);
    // lon2 = generate_random_double(-ngpt::DPI, ngpt::DPI);
    s = generate_random_double(0e0, ngpt::D2PI * R);
    az = generate_random_double(0e0, ngpt::D2PI);

    baz1 =
        ngpt::core::direct_vincenty(lat1, lon1, az, s, a, f, b, lat21, lon21);
  }
  auto t2 = std::chrono::high_resolution_clock::now();
  auto duration =
      std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
  std::cout << "\nImplementation 1: " << duration;

  t1 = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < 15000; ++i) {
    // random points 1 and 2
    lat1 = generate_random_double(-ngpt::DPI / 2e0, ngpt::DPI / 2e0);
    // lat2 = generate_random_double(-ngpt::DPI / 2e0, ngpt::DPI / 2e0);
    lon1 = generate_random_double(-ngpt::DPI, ngpt::DPI);
    // lon2 = generate_random_double(-ngpt::DPI, ngpt::DPI);
    s = generate_random_double(0e0, ngpt::D2PI * R);
    az = generate_random_double(0e0, ngpt::D2PI);

    baz2 =
        ngpt::core::direct_vincenty2(lat1, lon1, az, s, a, f, b, lat22, lon22);
  }
  t2 = std::chrono::high_resolution_clock::now();
  duration =
      std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
  std::cout << "\nImplementation 2: " << duration;

  // test from https://www.movable-type.co.uk/scripts/latlong-vincenty.html
  double p1_lat = hexd2rad(37, 57, 3.72030e0, -1);
  double p1_lon = hexd2rad(144, 25, 29.52440e0, 1);
  double bearing = hexd2rad(306, 52, 05.37e0, 1);
  double distance = 54972.271e0;
  baz1 = ngpt::core::direct_vincenty2(p1_lat, p1_lon, bearing, distance, a, f,
                                      b, lat21, lon21);
  baz2 = ngpt::core::direct_vincenty2(p1_lat, p1_lon, bearing, distance, a, f,
                                      b, lat22, lon22);
  int deg, min, sign;
  double sec;
  rad2hexd(p1_lat, deg, min, sec, sign);
  printf("\nStart Point: %03d %02d %9.5f", deg * sign, min, sec);
  rad2hexd(p1_lon, deg, min, sec, sign);
  printf(", %03d %02d %9.5f", deg * sign, min, sec);
  rad2hexd(bearing, deg, min, sec, sign);
  printf("\nBearing: %03d %02d %9.5f", deg * sign, min, sec);
  printf("\nDistance: %15.3f", distance);
  rad2hexd(lat21, deg, min, sec, sign);
  printf("\nDestination: %03d %02d %9.5f", deg * sign, min, sec);
  rad2hexd(lon21, deg, min, sec, sign);
  printf(", %03d %02d %9.5f", deg, min, sec);
  rad2hexd(baz1, deg, min, sec, sign);
  printf("\nFinal bearing: %03d %02d %9.5f", deg, min, sec);
  rad2hexd(lat22, deg, min, sec, sign);
  printf("\nDestination: %03d %02d %9.5f", deg * sign, min, sec);
  rad2hexd(lon22, deg, min, sec, sign);
  printf(", %03d %02d %9.5f", deg, min, sec);
  rad2hexd(baz2, deg, min, sec, sign);
  printf("\nFinal bearing: %03d %02d %9.5f", deg * sign, min, sec);

  printf("\n");
  return 0;
}
