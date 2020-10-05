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
double lat_at1, lat_at2, lat_at3;

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
    az2 = ngpt::core::direct_vincenty2(ar[0], ar[1], ar[2], ar[6], A, F, B,
                                       lat2, lon2);
    if (print)
      printf("\ndlat=%20.15f dlon=%20.15f dAz=%20.15f s=%15.3f",
             rad2seconds(ar[3] - lat2), rad2seconds(ar[4] - lon2),
             rad2seconds(ar[5] - az2), ar[6] * 1e-3);
    if (std::abs(rad2seconds(ar[3] - lat2)) > max_error_lat) {
      max_error_lat = std::abs(rad2seconds(ar[3] - lat2));
      lat_at1 = ar[3];
    }
    if (std::abs(rad2seconds(ar[4] - lon2)) > max_error_lon) {
      max_error_lon = std::abs(rad2seconds(ar[4] - lon2));
      lat_at2 = ar[3];
    }
    if (std::abs(rad2seconds(ar[5] - az2)) > max_error_az) {
      max_error_az = std::abs(rad2seconds(ar[5] - az2));
      lat_at3 = ar[3];
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

  for (const auto &batch : GeodTest_descr) {
    line_nr = 0;
    vec.clear();
    vec.reserve(batch.first);
    max_error_lat = std::numeric_limits<double>::min();
    max_error_lon = std::numeric_limits<double>::min();
    max_error_az = std::numeric_limits<double>::min();
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
    std::cout << "\n>> " << batch.second;
    std::cout
        << "\n-------------------------------------------------------------";
    printf("\nMax difs: lat:%20.15fsec lon:%20.15fsec az:%20.15fsec "
           "lat:%15.4f(m) lon:%15.4f(m)",
           max_error_lat, max_error_lon, max_error_az,
           infinitesimal_meridian_arc<ellipsoid::wgs84, double>(
               lat_at1, deg2rad(max_error_lat / 3600e0)),
           parallel_arc_length<ellipsoid::wgs84, double>(
               lat_at2, deg2rad(max_error_lon / 3600e0)));
  }
  printf("\nNumber of lines read: %15d", (int)line_count);

  printf("\n");
  return 0;
}
