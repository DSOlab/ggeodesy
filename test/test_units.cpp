#include "geodesy.hpp"
#include "test_help.hpp"
#include <cassert>
#include <iostream>
#include <random>

int main() {
  double angle_radians, angle_degrees, seconds, theta;
  int deg, min;

#ifdef CHECK_PRECISION
  double max_sin_diff = -std::numeric_limits<double>::max(),
    max_cos_diff = -std::numeric_limits<double>::max();
  double sin_diff, cos_diff;
#endif

  for (int i = 0; i < 150; ++i) {
    // radians to dec. degrees and back
    angle_radians = generate_random_double(-ngpt::D2PI, ngpt::D2PI);
    angle_degrees = ngpt::rad2deg(angle_radians);
    /* assert(std::abs(ngpt::deg2rad(angle_degrees) - angle_radians) <
     * ErrLimit); */
    assert(approxEqual(ngpt::deg2rad(angle_degrees), angle_radians));

    // dec. degrees to radians and back
    angle_degrees = generate_random_double(-360e0, 360e0);
    angle_radians = ngpt::deg2rad(angle_degrees);
    /* assert(std::abs(ngpt::rad2deg(angle_radians) - angle_degrees) <
     * ErrLimit); */
    assert(approxEqual(ngpt::rad2deg(angle_radians), angle_degrees));

    // dec. degrees to hexicondal degrees and back
    angle_degrees = generate_random_double(-360e0, 360e0);
    ngpt::decd2hexd(angle_degrees, deg, min, seconds);
    /*assert(std::abs(ngpt::hexd2decd(deg, min, seconds) - angle_degrees) <
           ErrLimit);*/
    assert(approxEqual(ngpt::hexd2decd(deg, min, seconds), angle_degrees));
    /*assert(std::abs(ngpt::hexd2rad(deg, min, seconds) -
                    ngpt::deg2rad(angle_degrees)) < ErrLimit);*/
    assert(approxEqual(ngpt::hexd2rad(deg, min, seconds),
                       ngpt::deg2rad(angle_degrees)));

    // radians to hexicondal degrees and back
    angle_radians = generate_random_double(-ngpt::D2PI, ngpt::D2PI);
    ngpt::rad2hexd(angle_radians, deg, min, seconds);
    /*assert(std::abs(ngpt::hexd2rad(deg, min, seconds) - angle_radians) <
           ErrLimit);*/
    assert(approxEqual(ngpt::hexd2rad(deg, min, seconds), angle_radians));
    /*assert(std::abs(ngpt::hexd2decd(deg, min, seconds) -
                    ngpt::rad2deg(angle_radians)) < ErrLimit);*/
    assert(approxEqual(ngpt::hexd2decd(deg, min, seconds),
                       ngpt::rad2deg(angle_radians)));

    // check the angle normalization
    angle_radians = generate_random_double(-ngpt::D2PI, ngpt::D2PI);
    theta = ngpt::normalize_angle(angle_radians, 0e0, ngpt::D2PI); // [0, 2π)
    /*assert(std::abs(std::sin(angle_radians) - std::sin(theta)) < ErrLimit);*/
#ifdef CHECK_PRECISION
    if ((sin_diff=std::abs(std::sin(angle_radians) - std::sin(theta)))>max_sin_diff)
      max_sin_diff = sin_diff;
    if ((cos_diff=std::abs(std::cos(angle_radians) - std::cos(theta)))>max_cos_diff)
      max_cos_diff = cos_diff;
#else
    assert(approxEqual(std::sin(angle_radians), std::sin(theta)));
    /*assert(std::abs(std::cos(angle_radians) - std::cos(theta)) < ErrLimit);*/
    assert(approxEqual(std::cos(angle_radians), std::cos(theta)));
#endif
    theta =
        ngpt::normalize_angle(angle_radians, -ngpt::DPI, ngpt::DPI); // [-π, π)
#ifdef CHECK_PRECISION
    if ((sin_diff=std::abs(std::sin(angle_radians) - std::sin(theta)))>max_sin_diff)
      max_sin_diff = sin_diff;
    if ((cos_diff=std::abs(std::cos(angle_radians) - std::cos(theta)))>max_cos_diff)
      max_cos_diff = cos_diff;
#else
    /*assert(std::abs(std::sin(angle_radians) - std::sin(theta)) < ErrLimit);*/
    assert(approxEqual(std::sin(angle_radians), std::sin(theta)));
    /*assert(std::abs(std::cos(angle_radians) - std::cos(theta)) < ErrLimit);*/
    assert(approxEqual(std::cos(angle_radians), std::cos(theta)));
#endif
  }

#ifdef CHECK_PRECISION
  printf("\nMax values for error:");
  printf("\n\tMax sinus diff   = %20.15f or %.15e",max_sin_diff, max_sin_diff);
  printf("\n\tMax cosinus diff = %20.15f or %.15e",max_cos_diff, max_cos_diff);
#endif

  return 0;
}
