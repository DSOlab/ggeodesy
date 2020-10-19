#include "geodesy.hpp"
#include "test_help.hpp"
#include "vincenty.hpp"
#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <fstream>
#include <iostream>
#include <random>

using namespace ngpt;

/// constants for wgs84 ellipsoid
constexpr double A = ellipsoid_traits<ellipsoid::wgs84>::a,
                 F = ellipsoid_traits<ellipsoid::wgs84>::f,
                 B = semi_minor<ellipsoid::wgs84>();

/// transform radians to seconds (of degrees)
double rad2seconds(double rad) noexcept {
  double ddeg = rad2deg(rad);
  return ddeg * 3600e0;
}

/// hold (absolute) max error values in seconds
double max_error_lat = std::numeric_limits<double>::min(),
       max_error_lon = std::numeric_limits<double>::min(),
       max_error_az = std::numeric_limits<double>::min();
double acc_error_lat, acc_error_lon, acc_error_az, acc_error_lat_m,
    acc_error_lon_m;
double lat_at1, lat_at2, lat_at3;
std::size_t failed_converges = 0;

/// struct to hold each line of the GeodTest file
/// see https://zenodo.org/record/32156
struct TestLine {
  // lat1[0], lon1[1], az1[2],
  // lat2[3], lon2[4], az2[5],
  // s12[6],  arc[7],  rs12[8],
  // t[9]
  double ar[10];

  void test_vincenty_direct(bool print = true) const {
    double lat2, lon2, az2;
    double err_lat, err_lon, err_az;
    az2 = ngpt::core::direct_vincenty2(ar[0], ar[1], ar[2], ar[6], A, F, B,
                                       lat2, lon2, 1e-15);
    if (print) {
      printf("\nInput:  lat=%+8.3f lon=%+8.3f Az=%+8.3f S=%8.1f",
             ngpt::rad2deg(ar[0]), ngpt::rad2deg(ar[1]), ngpt::rad2deg(ar[2]),
             ar[6] / 1000);
      printf("\nOutput: lon=%+8.3f should be %+8.3f", ngpt::rad2deg(lon2),
             ngpt::rad2deg(ar[4]));
    }
    if ((err_lat = std::abs(rad2seconds(ar[3] - lat2))) > max_error_lat) {
      max_error_lat = err_lat;
      lat_at1 = ar[3];
    }
    if ((err_lon = std::abs(rad2seconds(ar[4] - lon2))) > max_error_lon) {
      max_error_lon = err_lon;
      lat_at2 = ar[3];
    }
    if ((err_az = std::abs(rad2seconds(ar[5] - az2))) > max_error_az) {
      max_error_az = err_az;
      lat_at3 = ar[3];
    }
    acc_error_lat += err_lat;
    acc_error_lon += err_lon;
    acc_error_az += err_az;
    acc_error_lat_m += infinitesimal_meridian_arc<ellipsoid::wgs84, double>(
        lat2, deg2rad(err_lat / 3600e0));
    acc_error_lon_m += parallel_arc_length<ellipsoid::wgs84, double>(
        lat2, deg2rad(err_lon / 3600e0));
    return;
  }

  void test_vincenty_inverse(/*bool print = true*/) const {
    // err lat -> err in a12
    // err lon -> err in a21
    // err az  -> err in s12
    double s12, a12, a21;
    double err_lat, err_lon, err_az;
    try {
      s12 = ngpt::core::inverse_vincenty(ar[0], ar[1], ar[3], ar[4], A, F, B, a12,
          a21, 1e-12);
      if ((err_lat = std::abs(a12 - ar[2])) > max_error_lat) {
        max_error_lat = err_lat;
      }
      if ((err_lon = std::abs(a21 - ar[5])) > max_error_lon) {
        max_error_lon = err_lon;
      }
      if ((err_az = std::abs(s12 - ar[6])) > max_error_az) {
        max_error_az = err_az;
      }
      acc_error_lat += err_lat;
      acc_error_lon += err_lon;
      acc_error_az += err_az;
    } catch (std::out_of_range& e) {
      /*printf("\n[ERROR] Vincenty-Inverse failed to converge!");
      printf("\n        From point (%5.1f, %5.1f) to (%5.1f, %5.1f)",
          ngpt::rad2deg(ar[0]), ngpt::rad2deg(ar[1]), ngpt::rad2deg(ar[3]), 
          ngpt::rad2deg(ar[4]));*/
      ++failed_converges;
    }
    return;
  }
};

/// resolve a GeodTest line to a TestLine instance
int from_line(const char *line, TestLine &tl) {
  char *end;
  for (int i = 0; i < 10; i++) {
    tl.ar[i] = std::strtod(line, &end);
    assert(end != line);
    if (errno == ERANGE) {
      std::cerr << "\n[ERROR] Failed to decode line: \"" << line << "\"";
      return 1;
    }
    line = end;
  }
  for (int i = 0; i < 6; i++)
    tl.ar[i] = deg2rad(tl.ar[i]);
  tl.ar[7] = deg2rad(tl.ar[7]);
  return 0;
}

/// various constants
const char *GeodTest = "/home/xanthos/Builds/ggeodesy/test/GeodTest.dat";
constexpr std::size_t MAX_CHARS = 256;
using szstr_pair = std::pair<std::size_t, std::string>;
const std::vector<szstr_pair> GeodTest_descr = {
    {100000, "randomly distributed"},
    {50000, "nearly antipodal"},
    {50000, "short distances"},
    {50000, "one end near a pole"},
    {50000, "both ends near opposite poles"},
    {50000, "nearly meridional"},
    {50000, "nearly equatorial"},
    {50000, "running between vertices (a1 = a2 = 90deg)"},
    {50000, "ending close to vertices"}};
const std::size_t total_lines =
    std::accumulate(GeodTest_descr.begin(), GeodTest_descr.end(), 0,
                    [](const std::size_t &sum, const szstr_pair &pair) {
                      return sum + pair.first;
                    });
// assert(total_lines == 500000);

int main() {

  std::ifstream fin(GeodTest);
  if (!fin.is_open()) {
    std::cerr << "\n[ERROR] Failed to open file \"" << GeodTest << "\"";
    return 1;
  }

  std::size_t line_nr = 0, line_count = 0;
  char line[MAX_CHARS];
  std::vector<TestLine> vec;
  TestLine tl;

  /* testing direct vincenty */
  printf("\n%50s %17s %17s %17s %15s %15s", "Description", "Lat(seconds)",
         "Lon(seconds)", "Az(seconds)", "Lat(meters)", "Lon(meters)");
  std::cout
      << "\n-------------------------------------------------------------";
  for (const auto &batch : GeodTest_descr) {
    line_nr = 0;
    vec.clear();
    vec.reserve(batch.first);
    max_error_lat = std::numeric_limits<double>::min();
    max_error_lon = std::numeric_limits<double>::min();
    max_error_az = std::numeric_limits<double>::min();
    acc_error_lat = acc_error_lon = acc_error_az = acc_error_lat_m =
        acc_error_lon_m = 0e0;
    do {
      fin.getline(line, MAX_CHARS);
      if (!fin.good() || from_line(line, tl)) {
        return 1;
      }
      vec.emplace_back(tl);
    } while (++line_nr < batch.first);
    line_count += line_nr;
    for (const auto &ar : vec) {
      ar.test_vincenty_direct(false);
    }
    double mean_error_lat = acc_error_lat / (double)batch.first;
    double mean_error_lon = acc_error_lon / (double)batch.first;
    double mean_error_az = acc_error_az / (double)batch.first;
    double mean_error_lat_m = acc_error_lat_m / (double)batch.first;
    double mean_error_lon_m = acc_error_lon_m / (double)batch.first;
    printf("\n%50s %17.15f %17.15f %17.15f %15.4f %15.4f Max Err.",
           batch.second.c_str(), max_error_lat, max_error_lon, max_error_az,
           infinitesimal_meridian_arc<ellipsoid::wgs84, double>(
               lat_at1, deg2rad(max_error_lat / 3600e0)),
           parallel_arc_length<ellipsoid::wgs84, double>(
               lat_at2, deg2rad(max_error_lon / 3600e0)));
    printf("\n%50s %17.15f %17.15f %17.15f %15.4f %15.4f Mean Err.",
           batch.second.c_str(), mean_error_lat, mean_error_lon, mean_error_az,
           mean_error_lat_m, mean_error_lon_m);
  }
  printf("\nNumber of lines read: %15d", (int)line_count);

  /* testing inverse vincenty */
  fin.seekg(0);
  line_count = 0;
  printf("\n%50s %17s %17s %17s", "Description", "Az1->2(seconds)",
         "Az2->1(seconds)", "Distance(seconds)");
  std::cout
      << "\n-------------------------------------------------------------";
  for (const auto &batch : GeodTest_descr) {
    line_nr = 0;
    vec.clear();
    vec.reserve(batch.first);
    max_error_lat = std::numeric_limits<double>::min(); // a12
    max_error_lon = std::numeric_limits<double>::min(); // a21
    max_error_az = std::numeric_limits<double>::min();  // s12
    acc_error_lat = acc_error_lon = acc_error_az = 0e0;
    do {
      fin.getline(line, MAX_CHARS);
      if (!fin.good() || from_line(line, tl)) {
        return 1;
      }
      vec.emplace_back(tl);
    } while (++line_nr < batch.first);
    line_count += line_nr;
    bool doprint = false;
    for (const auto &ar : vec) {
      ar.test_vincenty_inverse();
    }
    double mean_error_lat = acc_error_lat / (double)batch.first;
    double mean_error_lon = acc_error_lon / (double)batch.first;
    double mean_error_az = acc_error_az / (double)batch.first;
    printf("\n%50s %17.15f %17.15f %17.15f Max Err.", batch.second.c_str(),
           max_error_lat, max_error_lon, max_error_az);
    printf("\n%50s %17.15f %17.15f %17.15f Mean Err.", batch.second.c_str(),
           mean_error_lat, mean_error_lon, mean_error_az);
  }
  printf("\nNumber of lines read: %15d failed to converge %3d%% (aka %15d)", 
      (int)line_count, (int)(failed_converges*100/line_count), (int)failed_converges);

  printf("\n");
  return 0;
}
