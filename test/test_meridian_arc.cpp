#include "ellipsoid.hpp"
#include "geodesy.hpp"
#include "meridian_arc_length.hpp"
#include "test_help.hpp"
#include <chrono>
#include <iostream>
#include <random>

/*
 * If CHECK_PRECISION compilation flag is set, then the program will check for
 * the precision between the meridian arc implementation (1, 2 and 3) and report
 * results. aka something like
 * Precision between implementation (1) and (2) is:        +0.0000323411m
 * Precision between implementation (1) and (3) is:        +0.0000012517m
 * Precision between implementation (2) and (3) is:        +0.0000311695m
 *
 * If the flag is not set, then the program will check (via assert) that the
 * results from the three different implementations agree to within a precision
 * of PRECISION meters.
 *
 * In either case, the program will also run a speed test on the three
 * implementations (aka measure and report execution time)
 */

#ifdef CHECK_PRECISION
double precision_12 = std::numeric_limits<double>::min();
double precision_13 = std::numeric_limits<double>::min();
double precision_23 = std::numeric_limits<double>::min();
#else
constexpr double PRECISION = 1e-4;
#endif

int main() {

  double angle_radians, lat;
  const double a = ngpt::ellipsoid_traits<ngpt::ellipsoid::grs80>::a;
  const double f = ngpt::ellipsoid_traits<ngpt::ellipsoid::grs80>::f;

  for (int i = 0; i < 10000; ++i) {
    lat = generate_random_double(-ngpt::DPI, ngpt::DPI);
    double s1 = ngpt::meridian_arc_length<ngpt::ellipsoid::grs80>(lat, 0);
    // printf("\nMeridian Arc Length to lat = %+5.1f is %+20.10fm",
    // ngpt::rad2deg(lat), s1);
    double s2 = ngpt::meridian_arc_length<ngpt::ellipsoid::grs80>(lat, 1);
    // printf("\nMeridian Arc Length to lat = %+5.1f is %+20.10fm",
    // ngpt::rad2deg(lat), s2);
    double s3 = ngpt::meridian_arc_length<ngpt::ellipsoid::grs80>(lat, 2);
    // printf("\nMeridian Arc Length to lat = %+5.1f is %+20.10fm",
    // ngpt::rad2deg(lat), s3);
#ifdef CHECK_PRECISION
    if (std::abs(s1 - s2) > precision_12)
      precision_12 = std::abs(s1 - s2);
    if (std::abs(s1 - s3) > precision_13)
      precision_13 = std::abs(s1 - s3);
    if (std::abs(s2 - s3) > precision_23)
      precision_23 = std::abs(s2 - s3);
#else
    assert(std::abs(s1 - s2) < PRECISION && std::abs(s2 - s3) < PRECISION);
#endif
  }

#ifdef CHECK_PRECISION
  printf("\nPrecision between implementation (1) and (2) is: %+20.10fm",
         precision_12);
  printf("\nPrecision between implementation (1) and (3) is: %+20.10fm",
         precision_13);
  printf("\nPrecision between implementation (2) and (3) is: %+20.10fm",
         precision_23);
#endif

  auto t1 = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < 10000; ++i) {
    lat = generate_random_double(-ngpt::DPI, ngpt::DPI);
    ngpt::core::meridian_arc_length_impl1(a, f, lat);
  }
  auto t2 = std::chrono::high_resolution_clock::now();
  auto duration =
      std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
  std::cout << "\nExecution time for impl1: " << duration << " microsec.";

  t1 = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < 10000; ++i) {
    lat = generate_random_double(-ngpt::DPI, ngpt::DPI);
    ngpt::core::meridian_arc_length_impl2(a, f, lat);
  }
  t2 = std::chrono::high_resolution_clock::now();
  duration =
      std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
  std::cout << "\nExecution time for impl2: " << duration << " microsec.";

  t1 = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < 10000; ++i) {
    lat = generate_random_double(-ngpt::DPI, ngpt::DPI);
    ngpt::core::meridian_arc_length_impl3(a, f, lat);
  }
  t2 = std::chrono::high_resolution_clock::now();
  duration =
      std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
  std::cout << "\nExecution time for impl3: " << duration << " microsec.";

  printf("\n");
  return 0;
}
