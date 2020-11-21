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
#ifdef DEBUG
#include <boost/geometry.hpp>
#include <boost/geometry/formulas/karney_direct.hpp>
#endif

using namespace ngpt;
using namespace boost::geometry;

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

/// hold (absolute) max error values in seconds
double mx_er_lat_v = std::numeric_limits<double>::min(), // max error for vincenty in rad
       mx_er_lon_v = std::numeric_limits<double>::min(),
       mx_er_az_v = std::numeric_limits<double>::min();
double mx_er_lat_k = std::numeric_limits<double>::min(), // max error for karney in rad
       mx_er_lon_k = std::numeric_limits<double>::min(),
       mx_er_az_k = std::numeric_limits<double>::min();
double mx_er_lat_bv = std::numeric_limits<double>::min(), // max error for boost
       mx_er_lon_bv = std::numeric_limits<double>::min(),
       mx_er_az_bv = std::numeric_limits<double>::min();
double mx_er_lat_bk = std::numeric_limits<double>::min(), // max error for boost
       mx_er_lon_bk = std::numeric_limits<double>::min(),
       mx_er_az_bk = std::numeric_limits<double>::min();
double ac_er_lat_v, ac_er_lon_v, ac_er_az_v,             // accumulated error for vincenty
  ac_er_lat_k, ac_er_lon_k, ac_er_az_k,                  // accumulated error for karney
  ac_er_lat_bv, ac_er_lon_bv, ac_er_az_bv,                  // accumulated error for boost
  ac_er_lat_bk, ac_er_lon_bk, ac_er_az_bk,                  // accumulated error for boost
  ac_er_lat_vm, ac_er_lon_vm,                            // accumulated error for vincenty in meters
  ac_er_lat_km, ac_er_lon_km,
  ac_er_lat_bvm, ac_er_lon_bvm,
  ac_er_lat_bkm, ac_er_lon_bkm;
double lat_at1, lat_at2, lat_at3;
// std::size_t failed_converges = 0;

/// struct to hold each line of the GeodTest file
/// see https://zenodo.org/record/32156
struct TestLine {
  // lat1[0], lon1[1], az1[2],
  // lat2[3], lon2[4], az2[5],
  // s12[6],  arc[7],  rs12[8],
  // t[9]
  double ar[10];

  void test_vincenty_direct() const {
    double lat2, lon2, az2;
    double err_lat, err_lon, err_az;
    az2 = ngpt::core::direct_vincenty2(ar[0], ar[1], ar[2], ar[6], A, F, B,
                                       lat2, lon2, 1e-15);
    if ((err_lat = std::abs(ar[3] - lat2)) > mx_er_lat_v) {
      mx_er_lat_v = err_lat;
      lat_at1 = ar[3];
    }
    if ((err_lon = std::abs(ar[4] - lon2)) > mx_er_lon_v) {
      mx_er_lon_v = err_lon;
      lat_at2 = ar[3];
    }
    if ((err_az = std::abs(ar[5] - az2)) > mx_er_az_v) {
      mx_er_az_v = err_az;
      lat_at3 = ar[3];
    }
    ac_er_lat_v += err_lat;
    ac_er_lon_v += err_lon;
    ac_er_az_v += err_az;
    ac_er_lat_vm += infinitesimal_meridian_arc<ellipsoid::wgs84, double>(
        lat2, deg2rad(err_lat / 3600e0));
    ac_er_lon_vm += parallel_arc_length<ellipsoid::wgs84, double>(
        lat2, deg2rad(err_lon / 3600e0));
  }
  void test_boostv_direct() const {
    double lat2, lon2, az2;
    double err_lat, err_lon, err_az;
    result_drc =
        vincenty_direct_type::apply(ar[1], ar[0], ar[6], ar[2], spheroid);
    if ((err_lat = std::abs(ar[3] - result_drc.lat2)) > mx_er_lat_bv) {
      mx_er_lat_bv = err_lat;
    }
    if ((err_lon = std::abs(ar[4] - result_drc.lon2)) >
        mx_er_lon_bv) {
      mx_er_lon_bv = err_lon;
    }
    if ((err_az = std::abs(ar[5] - result_drc.reverse_azimuth)) > mx_er_az_bv) {
      mx_er_az_bv = err_az;
    }
    ac_er_lat_bv += err_lat;
    ac_er_lon_bv += err_lon;
    ac_er_az_bv += err_az;
    ac_er_lat_bvm += infinitesimal_meridian_arc<ellipsoid::wgs84, double>(
        result_drc.lat2, deg2rad(err_lat / 3600e0));
    ac_er_lon_bvm += parallel_arc_length<ellipsoid::wgs84, double>(
        result_drc.lat2, deg2rad(err_lon / 3600e0));
  }
  void test_boostk_direct() const {
    double lat2, lon2, az2;
    double err_lat, err_lon, err_az;
    result_drc =
        karney_direct_type::apply(ar[1], ar[0], ar[6], ar[2], spheroid);
    if ((err_lat = std::abs(ar[3] - result_drc.lat2)) > mx_er_lat_bk) {
      mx_er_lat_bk = err_lat;
    }
    if ((err_lon = std::abs(ar[4] - result_drc.lon2)) >
        mx_er_lon_bk) {
      mx_er_lon_bk = err_lon;
    }
    if ((err_az = std::abs(ar[5] - result_drc.reverse_azimuth)) > mx_er_az_bk) {
      mx_er_az_bk = err_az;
    }
    ac_er_lat_bk += err_lat;
    ac_er_lon_bk += err_lon;
    ac_er_az_bk += err_az;
    ac_er_lat_bkm += infinitesimal_meridian_arc<ellipsoid::wgs84, double>(
        result_drc.lat2, deg2rad(err_lat / 3600e0));
    ac_er_lon_bkm += parallel_arc_length<ellipsoid::wgs84, double>(
        result_drc.lat2, deg2rad(err_lon / 3600e0));
  }
  void test_karney_direct() const {
    double lat2, lon2, az2;
    double err_lat, err_lon, err_az;
    // check against karney
    az2 = ngpt::core::direct_karney(ar[0], ar[1], ar[2], ar[6], A, F, B,
                                       lat2, lon2);
    if ((err_lat = std::abs(ar[3] - lat2)) > mx_er_lat_k) {
      mx_er_lat_k = err_lat;
      lat_at1 = ar[3];
    }
    if ((err_lon = std::abs(ar[4] - lon2)) > mx_er_lon_k) {
      mx_er_lon_k = err_lon;
      lat_at2 = ar[3];
    }
    if ((err_az = std::abs(ar[5] - az2)) > mx_er_az_k) {
      mx_er_az_k = err_az;
      lat_at3 = ar[3];
    }
    ac_er_lat_k += err_lat;
    ac_er_lon_k += err_lon;
    ac_er_az_k += err_az;
    ac_er_lat_km += infinitesimal_meridian_arc<ellipsoid::wgs84, double>(
        lat2, deg2rad(err_lat / 3600e0));
    ac_er_lon_km += parallel_arc_length<ellipsoid::wgs84, double>(
        lat2, deg2rad(err_lon / 3600e0));
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
    mx_er_lat_v = std::numeric_limits<double>::min(), // max error for vincenty in rad
                mx_er_lon_v = std::numeric_limits<double>::min(),
                mx_er_az_v = std::numeric_limits<double>::min();
    mx_er_lat_k = std::numeric_limits<double>::min(), // max error for karney in rad
                mx_er_lon_k = std::numeric_limits<double>::min(),
                mx_er_az_k = std::numeric_limits<double>::min();
    mx_er_lat_bv = std::numeric_limits<double>::min(), // max error for boost
                mx_er_lon_bv = std::numeric_limits<double>::min(),
                mx_er_az_bv = std::numeric_limits<double>::min();
    mx_er_lat_bk = std::numeric_limits<double>::min(), // max error for boost
                mx_er_lon_bk = std::numeric_limits<double>::min(),
                mx_er_az_bk = std::numeric_limits<double>::min();
    ac_er_lat_v= ac_er_lon_v= ac_er_az_v=             // accumulated error for vincenty
      ac_er_lat_k= ac_er_lon_k= ac_er_az_k=                  // accumulated error for karney
      ac_er_lat_bv= ac_er_lon_bv= ac_er_az_bv=                  // accumulated error for boost
      ac_er_lat_bk= ac_er_lon_bk= ac_er_az_bk=                  // accumulated error for boost
      ac_er_lat_vm= ac_er_lon_vm=                            // accumulated error for vincenty in meters
      ac_er_lat_km= ac_er_lon_km=
      ac_er_lat_bvm= ac_er_lon_bvm=
      ac_er_lat_bkm= ac_er_lon_bkm = 0e0;
    do {
      fin.getline(line, MAX_CHARS);
      if (!fin.good() || from_line(line, tl)) {
        return 1;
      }
      vec.emplace_back(tl);
    } while (++line_nr < batch.first);
    line_count += line_nr;
    for (const auto &ar : vec) {
      ar.test_vincenty_direct();
      ar.test_karney_direct();
      ar.test_boostv_direct();
      ar.test_boostk_direct();
    }
    double mean_error_lat = ac_er_lat_v / (double)batch.first;
    double mean_error_lon = ac_er_lon_v / (double)batch.first;
    double mean_error_az = ac_er_az_v / (double)batch.first;
    double mean_error_lat_m = ac_er_lat_vm / (double)batch.first;
    double mean_error_lon_m = ac_er_lon_vm / (double)batch.first;
    printf("\n%50s %17.15f %17.15f %17.15f %15.4f %15.4f Max Err. Vincenty",
           batch.second.c_str(), rad2sec(mx_er_lat_v), rad2sec(mx_er_lon_v), rad2sec(mx_er_az_v),
           infinitesimal_meridian_arc<ellipsoid::wgs84, double>(
               lat_at1, mx_er_lat_v),
           parallel_arc_length<ellipsoid::wgs84, double>(
               lat_at2, mx_er_lon_v));
    printf("\n%50s %17.15f %17.15f %17.15f %15.4f %15.4f Mean Err. Vincenty",
           batch.second.c_str(), rad2sec(mean_error_lat), rad2sec(mean_error_lon), rad2sec(mean_error_az),
           mean_error_lat_m, mean_error_lon_m);
    mean_error_lat = ac_er_lat_k / (double)batch.first;
    mean_error_lon = ac_er_lon_k / (double)batch.first;
    mean_error_az = ac_er_az_k / (double)batch.first;
    mean_error_lat_m = ac_er_lat_km / (double)batch.first;
    mean_error_lon_m = ac_er_lon_km / (double)batch.first;
    printf("\n%50s %17.15f %17.15f %17.15f %15.4f %15.4f Max Err. Karney",
           batch.second.c_str(), rad2sec(mx_er_lat_k), rad2sec(mx_er_lon_k), rad2sec(mx_er_az_k),
           infinitesimal_meridian_arc<ellipsoid::wgs84, double>(
               lat_at1, mx_er_lat_k),
           parallel_arc_length<ellipsoid::wgs84, double>(
               lat_at2, mx_er_lon_k));
    printf("\n%50s %17.15f %17.15f %17.15f %15.4f %15.4f Mean Err. Karney",
           batch.second.c_str(), rad2sec(mean_error_lat), rad2sec(mean_error_lon), rad2sec(mean_error_az),
           mean_error_lat_m, mean_error_lon_m);
    mean_error_lat = ac_er_lat_bv / (double)batch.first;
    mean_error_lon = ac_er_lon_bv / (double)batch.first;
    mean_error_az = ac_er_az_bv / (double)batch.first;
    mean_error_lat_m = ac_er_lat_bvm / (double)batch.first;
    mean_error_lon_m = ac_er_lon_bvm / (double)batch.first;
    printf("\n%50s %17.15f %17.15f %17.15f %15.4f %15.4f Max Err. Boost (Vincenty)",
           batch.second.c_str(), rad2sec(mx_er_lat_bv), rad2sec(mx_er_lon_bv), rad2sec(mx_er_az_bv),
           infinitesimal_meridian_arc<ellipsoid::wgs84, double>(
               lat_at1, mx_er_lat_bv),
           parallel_arc_length<ellipsoid::wgs84, double>(
               lat_at2, mx_er_lon_bv));
    printf("\n%50s %17.15f %17.15f %17.15f %15.4f %15.4f Mean Err. Boost (Vincenty)",
           batch.second.c_str(), rad2sec(mean_error_lat), rad2sec(mean_error_lon), rad2sec(mean_error_az),
           mean_error_lat_m, mean_error_lon_m);
    mean_error_lat = ac_er_lat_bk / (double)batch.first;
    mean_error_lon = ac_er_lon_bk / (double)batch.first;
    mean_error_az = ac_er_az_bk / (double)batch.first;
    mean_error_lat_m = ac_er_lat_bkm / (double)batch.first;
    mean_error_lon_m = ac_er_lon_bkm / (double)batch.first;
    printf("\n%50s %17.15f %17.15f %17.15f %15.4f %15.4f Max Err. Boost (Karney)",
           batch.second.c_str(), rad2sec(mx_er_lat_bk), rad2sec(mx_er_lon_bk), rad2sec(mx_er_az_bk),
           infinitesimal_meridian_arc<ellipsoid::wgs84, double>(
               lat_at1, mx_er_lat_bk),
           parallel_arc_length<ellipsoid::wgs84, double>(
               lat_at2, mx_er_lon_bk));
    printf("\n%50s %17.15f %17.15f %17.15f %15.4f %15.4f Mean Err. Boost (Karney)",
           batch.second.c_str(), rad2sec(mean_error_lat), rad2sec(mean_error_lon), rad2sec(mean_error_az),
           mean_error_lat_m, mean_error_lon_m);
  }
  printf("\nNumber of lines read: %15d", (int)line_count);

  printf("\n");
  return 0;
}
