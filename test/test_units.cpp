#include "geodesy.hpp"
#include "test_help.hpp"
#include <cassert>
#include <iostream>
#include <random>

int main() {
  double angle_radians, angle_degrees, seconds, theta;
  int deg, min, sign;

#ifdef CHECK_PRECISION
  double max_sin_diff = -std::numeric_limits<double>::max(),
         max_cos_diff = -std::numeric_limits<double>::max();
  double sin_diff, cos_diff;
  double hex2decd_max = -std::numeric_limits<double>::max(),
         hexd2rad_max = -std::numeric_limits<double>::max();
  double hex2decd_max_angle = 0e0, hexd2rad_max_angle = 0e0;
#else
  constexpr double PRECISION = 1e-15;
#endif

  // test hexicondal degrees sign
  int deg1, min1, deg2, min2, sgn1, sgn2;
  double sec1, sec2, a1, a2;
  ngpt::decd2hexd(10e0, deg1, min1, sec1, sgn1);
  ngpt::decd2hexd(-10e0, deg2, min2, sec2, sgn2);
  assert(deg1 == deg2 && (min1 == min2 && sec1 == sec2));
  assert(sgn1 == -sgn2);
  //
  a1 = ngpt::hexd2decd(deg1, min1, sec1, sgn1);
  a2 = ngpt::hexd2decd(deg2, min2, sec2, sgn2);
  assert(a1 == -a2);
  // note the sign of deg parameter is not considered!
  a2 = ngpt::hexd2decd(-deg1, min1, sec1, sgn1);
  assert(a1 == a2);
  //
  a1 = ngpt::hexd2rad(deg1, min1, sec1, sgn1);
  a2 = ngpt::hexd2rad(deg2, min2, sec2, sgn2);
  assert(a1 == -a2);
  // note the sign of deg parameter is not considered!
  a1 = ngpt::hexd2rad(-deg1, min1, sec1, sgn1);
  assert(a1 == -a2);

  for (int i = 0; i < 1500; ++i) {
    // radians to dec. degrees and back
    angle_radians = generate_random_double(-ngpt::D2PI, ngpt::D2PI);
    angle_degrees = ngpt::rad2deg(angle_radians);
    assert(approxEqual(ngpt::deg2rad(angle_degrees), angle_radians));

    // dec. degrees to radians and back
    angle_degrees = generate_random_double(-360e0, 360e0);
    angle_radians = ngpt::deg2rad(angle_degrees);
    assert(approxEqual(ngpt::rad2deg(angle_radians), angle_degrees));

    // dec. degrees to hexicondal degrees and back
    angle_degrees = generate_random_double(-360e0, 360e0);
    ngpt::decd2hexd(angle_degrees, deg, min, seconds, sign);
#ifdef CHECK_PRECISION
    if (!approxEqual(ngpt::hexd2decd(deg, min, seconds, sign), angle_degrees)) {
      double diff =
          std::abs(ngpt::hexd2decd(deg, min, seconds, sign) - angle_degrees);
      if (hex2decd_max < diff) {
        printf("\nhexd2decd: degrees: %+20.15f, transformation: %+20.15f",
               angle_degrees, ngpt::hexd2decd(deg, min, seconds, sign));
        hex2decd_max = diff;
        hex2decd_max_angle = angle_degrees;
      }
    }
#else
    assert(
        approxEqual(ngpt::hexd2decd(deg, min, seconds, sign), angle_degrees));
    assert(approxEqual(ngpt::hexd2rad(deg, min, seconds, sign),
                       ngpt::deg2rad(angle_degrees)));
#endif

    // radians to hexicondal degrees and back
    angle_radians = generate_random_double(-ngpt::D2PI, ngpt::D2PI);
    ngpt::rad2hexd(angle_radians, deg, min, seconds, sign);
#ifdef CHECK_PRECISION
    if (!approxEqual(ngpt::hexd2rad(deg, min, seconds, sign), angle_radians)) {
      if (hexd2rad_max <
          std::abs(ngpt::hexd2rad(deg, min, seconds, sign) - angle_radians)) {
        hexd2rad_max =
            std::abs(ngpt::hexd2rad(deg, min, seconds, sign) - angle_radians);
        hexd2rad_max_angle = ngpt::rad2deg(angle_radians);
      }
    }
#else
    assert(approxEqual(ngpt::hexd2rad(deg, min, seconds, sign), angle_radians));
    assert(approxEqual(ngpt::hexd2decd(deg, min, seconds, sign),
                       ngpt::rad2deg(angle_radians)));
#endif

    // check the angle normalization
    angle_radians = generate_random_double(-ngpt::D2PI, ngpt::D2PI);
    theta = ngpt::normalize_angle(angle_radians, 0e0, ngpt::D2PI); // [0, 2π)
#ifdef CHECK_PRECISION
    if ((sin_diff = std::abs(std::sin(angle_radians) - std::sin(theta))) >
        max_sin_diff)
      max_sin_diff = sin_diff;
    if ((cos_diff = std::abs(std::cos(angle_radians) - std::cos(theta))) >
        max_cos_diff)
      max_cos_diff = cos_diff;
#else
    /* the following fails on some rare occasions
     * assert(approxEqual(std::sin(angle_radians), std::sin(theta)));
     * assert(approxEqual(std::cos(angle_radians), std::cos(theta)));
     */
    assert(std::abs(std::sin(angle_radians) - std::sin(theta)) < PRECISION);
    assert(std::abs(std::cos(angle_radians) - std::cos(theta)) < PRECISION);
#endif
    theta =
        ngpt::normalize_angle(angle_radians, -ngpt::DPI, ngpt::DPI); // [-π, π)
#ifdef CHECK_PRECISION
    if ((sin_diff = std::abs(std::sin(angle_radians) - std::sin(theta))) >
        max_sin_diff)
      max_sin_diff = sin_diff;
    if ((cos_diff = std::abs(std::cos(angle_radians) - std::cos(theta))) >
        max_cos_diff)
      max_cos_diff = cos_diff;
#else
    /* the following may fail on rare occasions
     * assert(approxEqual(std::sin(angle_radians), std::sin(theta)));
     * assert(approxEqual(std::cos(angle_radians), std::cos(theta)));
     */
    assert(std::abs(std::sin(angle_radians) - std::sin(theta)) < PRECISION);
    assert(std::abs(std::cos(angle_radians) - std::cos(theta)) < PRECISION);
#endif
  }

#ifdef CHECK_PRECISION
  printf("\nMax values for error:");
  printf("\n\tMax hex2decd = %20.15f deg. or %.15e at angle=%.15e degrees",
         hex2decd_max, hex2decd_max, hex2decd_max_angle);
  printf("\n\tMax hex2rad  = %20.15f rad  or %.15e at angle=%.15e degrees",
         hexd2rad_max, hexd2rad_max, hexd2rad_max_angle);
  printf("\n\tMax sinus diff   = %20.15f or %.15e", max_sin_diff, max_sin_diff);
  printf("\n\tMax cosinus diff = %20.15f or %.15e", max_cos_diff, max_cos_diff);
#endif

  printf("\n");
  return 0;
}
