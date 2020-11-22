#include "geodesy.hpp"
#include "vincenty.hpp"
#include "units.hpp"
#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <fstream>
#include <iostream>
#include <random>
#include <boost/geometry.hpp>
#include <boost/geometry/formulas/karney_direct.hpp>
#include <chrono>

using namespace ngpt;
using namespace boost::geometry;
using namespace std::chrono;

/// constants for wgs84 ellipsoid
constexpr double A = ellipsoid_traits<ellipsoid::wgs84>::a,
                 F = ellipsoid_traits<ellipsoid::wgs84>::f,
                 B = semi_minor<ellipsoid::wgs84>();

// For storing the resulting values.
formula::result_direct<double> result_drc;

// WGS-84 spheroid.
srs::spheroid<double> spheroid(A, B);

// Define the strategy.
typedef formula::vincenty_direct<double, true, true, false, false>
    vincenty_direct_type;
typedef formula::karney_direct<double, true, true, false, false, 6u>
    karney_direct_type;

/// struct to hold each line of the GeodTest file
/// see https://zenodo.org/record/32156
struct TestLine {
  // lat1[0], lon1[1], az1[2],
  // lat2[3], lon2[4], az2[5],
  // s12[6],  arc[7],  rs12[8],
  // t[9]
  double ar[10];
  double lat1() const noexcept { return ar[0]; }
  double lon1() const noexcept { return ar[1]; }
  double lat2() const noexcept { return ar[3]; }
  double lon2() const noexcept { return ar[4]; }
  double forward_azimuth() const noexcept { return ar[2]; }
  double backward_azimuth() const noexcept { return ar[5]; }
  double distance() const noexcept { return ar[6]; }
}; // TestLine

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
// const char *GeodTest = "/home/xanthos/Software/ggeodesy/test/GeodTest.dat";
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

auto test_ggeodesy_vincenty(const std::vector<TestLine> &vec, double* maxer, double* meaner) {
  double maxer_lat=std::numeric_limits<double>::min(),
    maxer_lon=std::numeric_limits<double>::min(),
    maxer_raz=std::numeric_limits<double>::min();
  double meaner_lat=0e0, meaner_lon=0e0, meaner_raz=0e0;
  double az2, lat2, lon2, err_lat, err_lon, err_az;
  auto start = high_resolution_clock::now();
  for (const auto& dataset : vec) {
    az2 = ngpt::core::direct_vincenty2(dataset.lat1(), dataset.lon1(), dataset.forward_azimuth(), dataset.distance(), A, F, B,
                                       lat2, lon2, 1e-15);
    if ((err_lat = std::abs(dataset.lat2() - lat2)) > maxer_lat) {
      maxer_lat = err_lat;
    }
    if ((err_lon = std::abs(dataset.lon2() - lon2)) > maxer_lon) {
      maxer_lon = err_lon;
    }
    if ((err_az = std::abs(dataset.backward_azimuth() - az2)) > maxer_raz) {
      maxer_raz = err_az;
    }
    meaner_lat += err_lat;
    meaner_lon += err_lon;
    meaner_raz += err_az;
  }
  auto stop = high_resolution_clock::now();
  maxer[0] = maxer_lat;
  maxer[1] = maxer_lon;
  maxer[2] = maxer_raz;
  meaner[0] = meaner_lat / vec.size();
  meaner[1] = meaner_lon / vec.size();
  meaner[2] = meaner_raz / vec.size();
  return duration_cast<microseconds>(stop - start);
}
auto test_ggeodesy_karney(const std::vector<TestLine> &vec, double* maxer, double* meaner) {
  double maxer_lat=std::numeric_limits<double>::min(),
    maxer_lon=std::numeric_limits<double>::min(),
    maxer_raz=std::numeric_limits<double>::min();
  double meaner_lat=0e0, meaner_lon=0e0, meaner_raz=0e0;
  double az2, lat2, lon2, err_lat, err_lon, err_az;
  auto start = high_resolution_clock::now();
  for (const auto& dataset : vec) {
    az2 = ngpt::core::direct_karney(dataset.lat1(), dataset.lon1(), dataset.forward_azimuth(), dataset.distance(), A, F, B,
                                       lat2, lon2);
    if ((err_lat = std::abs(dataset.lat2() - lat2)) > maxer_lat) {
      maxer_lat = err_lat;
    }
    if ((err_lon = std::abs(dataset.lon2() - lon2)) > maxer_lon) {
      maxer_lon = err_lon;
    }
    if ((err_az = std::abs(dataset.backward_azimuth() - az2)) > maxer_raz) {
      maxer_raz = err_az;
    }
    meaner_lat += err_lat;
    meaner_lon += err_lon;
    meaner_raz += err_az;
  }
  auto stop = high_resolution_clock::now();
  maxer[0] = maxer_lat;
  maxer[1] = maxer_lon;
  maxer[2] = maxer_raz;
  meaner[0] = meaner_lat / vec.size();
  meaner[1] = meaner_lon / vec.size();
  meaner[2] = meaner_raz / vec.size();
  return duration_cast<microseconds>(stop - start);
}
auto test_boost_vincenty(const std::vector<TestLine> &vec, double* maxer, double* meaner) {
  double maxer_lat=std::numeric_limits<double>::min(),
    maxer_lon=std::numeric_limits<double>::min(),
    maxer_raz=std::numeric_limits<double>::min();
  double meaner_lat=0e0, meaner_lon=0e0, meaner_raz=0e0;
  double az2, lat2, lon2, err_lat, err_lon, err_az;
  auto start = high_resolution_clock::now();
  for (const auto& dataset : vec) {
    result_drc =
        vincenty_direct_type::apply(dataset.lon1(), dataset.lat1(), dataset.distance(), dataset.forward_azimuth(), spheroid);
    if ((err_lat = std::abs(dataset.lat2() - result_drc.lat2)) > maxer_lat) {
      maxer_lat = err_lat;
    }
    if ((err_lon = std::abs(dataset.lon2() - result_drc.lon2)) > maxer_lon) {
      maxer_lon = err_lon;
    }
    if ((err_az = std::abs(dataset.backward_azimuth() - result_drc.reverse_azimuth)) > maxer_raz) {
      maxer_raz = err_az;
    }
    meaner_lat += err_lat;
    meaner_lon += err_lon;
    meaner_raz += err_az;
  }
  auto stop = high_resolution_clock::now();
  maxer[0] = maxer_lat;
  maxer[1] = maxer_lon;
  maxer[2] = maxer_raz;
  meaner[0] = meaner_lat / vec.size();
  meaner[1] = meaner_lon / vec.size();
  meaner[2] = meaner_raz / vec.size();
  return duration_cast<microseconds>(stop - start);
}
auto test_boost_karney(const std::vector<TestLine> &cvec, double* maxer, double* meaner) {
  double maxer_lat=std::numeric_limits<double>::min(),
    maxer_lon=std::numeric_limits<double>::min(),
    maxer_raz=std::numeric_limits<double>::min();
  double meaner_lat=0e0, meaner_lon=0e0, meaner_raz=0e0;
  double az2, lat2, lon2, err_lat, err_lon, err_az;
  auto vec = cvec;
  for (auto& dataset : vec) {
    for (int i=0; i<6; i++) dataset.ar[i] = rad2deg(dataset.ar[i]);
  }
  auto start = high_resolution_clock::now();
  for (const auto& dataset : vec) {
    result_drc =
        karney_direct_type::apply(dataset.lon1(), dataset.lat1(), dataset.distance(), dataset.forward_azimuth(), spheroid);
    if ((err_lat = std::abs(dataset.lat2() - result_drc.lat2)) > maxer_lat) {
      maxer_lat = err_lat;
    }
    if ((err_lon = std::abs(dataset.lon2() - result_drc.lon2)) > maxer_lon) {
      maxer_lon = err_lon;
    }
    if ((err_az = std::abs(dataset.backward_azimuth() - result_drc.reverse_azimuth)) > maxer_raz) {
      maxer_raz = err_az;
    }
    meaner_lat += err_lat;
    meaner_lon += err_lon;
    meaner_raz += err_az;
  }
  auto stop = high_resolution_clock::now();
  maxer[0] = deg2rad(maxer_lat);
  maxer[1] = deg2rad(maxer_lon);
  maxer[2] = deg2rad(maxer_raz);
  meaner[0] = deg2rad(meaner_lat) / vec.size();
  meaner[1] = deg2rad(meaner_lon) / vec.size();
  meaner[2] = deg2rad(meaner_raz) / vec.size();
  return duration_cast<microseconds>(stop - start);
}

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
  double maxer[3], meaner[3];

  /* */
  printf("\n%50s %17s %17s %17s %15s %15s", "Description", "Lat(seconds)",
         "Lon(seconds)", "Az(seconds)", "Lat(meters)", "Lon(meters)");
  std::cout
      << "\n-------------------------------------------------------------";
  printf("> ggeodesy Vincenty");
  for (const auto &batch : GeodTest_descr) {
    line_nr = 0;
    vec.clear();
    vec.reserve(batch.first);
    std::size_t running_time = 0;
    do {
      fin.getline(line, MAX_CHARS);
      if (!fin.good() || from_line(line, tl)) {
        return 1;
      }
      vec.emplace_back(tl);
    } while (++line_nr < batch.first);
    line_count += line_nr;
    auto duration = test_ggeodesy_vincenty(vec, maxer, meaner);
    printf("\n%50s %17.15f %17.15f %17.15f Max Err.",
           batch.second.c_str(), rad2sec(maxer[0]), rad2sec(maxer[1]), rad2sec(maxer[2]));
    printf("\n%50s %17.15f %17.15f %17.15f Mean Err. (time: %10ld microsec)",
           batch.second.c_str(), rad2sec(meaner[0]), rad2sec(meaner[1]), rad2sec(meaner[2]), duration.count());
    running_time += duration.count();
  }
  fin.seekg(0);
  std::cout
      << "\n-------------------------------------------------------------";
  printf("> ggeodesy Karney");
  for (const auto &batch : GeodTest_descr) {
    line_nr = 0;
    vec.clear();
    vec.reserve(batch.first);
    std::size_t running_time = 0;
    do {
      fin.getline(line, MAX_CHARS);
      if (!fin.good() || from_line(line, tl)) {
        return 1;
      }
      vec.emplace_back(tl);
    } while (++line_nr < batch.first);
    line_count += line_nr;
    auto duration = test_ggeodesy_karney(vec, maxer, meaner);
    printf("\n%50s %17.15f %17.15f %17.15f Max Err.",
           batch.second.c_str(), rad2sec(maxer[0]), rad2sec(maxer[1]), rad2sec(maxer[2]));
    printf("\n%50s %17.15f %17.15f %17.15f Mean Err. (time: %10ld microsec)",
           batch.second.c_str(), rad2sec(meaner[0]), rad2sec(meaner[1]), rad2sec(meaner[2]), duration.count());
    running_time += duration.count();
  }
  fin.seekg(0);
  std::cout
      << "\n-------------------------------------------------------------";
  printf("> Boost Vincenty");
  for (const auto &batch : GeodTest_descr) {
    line_nr = 0;
    vec.clear();
    vec.reserve(batch.first);
    std::size_t running_time = 0;
    do {
      fin.getline(line, MAX_CHARS);
      if (!fin.good() || from_line(line, tl)) {
        return 1;
      }
      vec.emplace_back(tl);
    } while (++line_nr < batch.first);
    line_count += line_nr;
    auto duration = test_boost_vincenty(vec, maxer, meaner);
    printf("\n%50s %17.15f %17.15f %17.15f Max Err.",
           batch.second.c_str(), rad2sec(maxer[0]), rad2sec(maxer[1]), rad2sec(maxer[2]));
    printf("\n%50s %17.15f %17.15f %17.15f Mean Err. (time: %10ld microsec)",
           batch.second.c_str(), rad2sec(meaner[0]), rad2sec(meaner[1]), rad2sec(meaner[2]), duration.count());
    running_time += duration.count();
  }
  fin.seekg(0);
  std::cout
      << "\n-------------------------------------------------------------";
  printf("> Boost Karney");
  for (const auto &batch : GeodTest_descr) {
    line_nr = 0;
    vec.clear();
    vec.reserve(batch.first);
    std::size_t running_time = 0;
    do {
      fin.getline(line, MAX_CHARS);
      if (!fin.good() || from_line(line, tl)) {
        return 1;
      }
      vec.emplace_back(tl);
    } while (++line_nr < batch.first);
    line_count += line_nr;
    auto duration = test_boost_karney(vec, maxer, meaner);
    printf("\n%50s %17.15f %17.15f %17.15f Max Err.",
           batch.second.c_str(), rad2sec(maxer[0]), rad2sec(maxer[1]), rad2sec(maxer[2]));
    printf("\n%50s %17.15f %17.15f %17.15f Mean Err. (time: %10ld microsec)",
           batch.second.c_str(), rad2sec(meaner[0]), rad2sec(meaner[1]), rad2sec(meaner[2]), duration.count());
    running_time += duration.count();
  }
  printf("\nNumber of lines read: %15d", (int)line_count);

  printf("\n");
  return 0;
}
