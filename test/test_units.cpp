#include "geodesy.hpp"
#include "test_help.hpp"
#include "units.hpp"
#include <cassert>
#include <iostream>
#include <random>

int main() {
  double angle_radians, angle_degrees, seconds, theta;
  int deg, min, sign;

  // test hexicondal degrees sign
  int deg1, min1, deg2, min2, sgn1, sgn2;
  double sec1, sec2, a1, a2;
  dso::decd2hexd(10e0, deg1, min1, sec1, sgn1);
  dso::decd2hexd(-10e0, deg2, min2, sec2, sgn2);
  assert(deg1 == deg2 && (min1 == min2 && sec1 == sec2));
  assert(sgn1 == -sgn2);
  //
  a1 = dso::hexd2decd(deg1, min1, sec1, sgn1);
  a2 = dso::hexd2decd(deg2, min2, sec2, sgn2);
  assert(a1 == -a2);
  // note the sign of deg parameter is not considered!
  a2 = dso::hexd2decd(-deg1, min1, sec1, sgn1);
  assert(a1 == a2);
  //
  a1 = dso::hexd2rad(deg1, min1, sec1, sgn1);
  a2 = dso::hexd2rad(deg2, min2, sec2, sgn2);
  assert(a1 == -a2);
  // note the sign of deg parameter is not considered!
  a1 = dso::hexd2rad(-deg1, min1, sec1, sgn1);
  assert(a1 == -a2);

  for (int i = 0; i < 15000; ++i) {
    // radians to dec. degrees and back
    angle_radians = generate_random_double(-dso::D2PI, dso::D2PI);
    angle_degrees = dso::rad2deg(angle_radians);
    assert(approxEqual(dso::deg2rad(angle_degrees), angle_radians));

    // dec. degrees to radians and back
    angle_degrees = generate_random_double(-360e0, 360e0);
    angle_radians = dso::deg2rad(angle_degrees);
    assert(approxEqual(dso::rad2deg(angle_radians), angle_degrees));

    // dec. degrees to hexicondal degrees and back
    angle_degrees = generate_random_double(-360e0, 360e0);
    dso::decd2hexd(angle_degrees, deg, min, seconds, sign);
    assert(
        approxEqual(dso::hexd2decd(deg, min, seconds, sign), angle_degrees));
    assert(approxEqual(dso::hexd2rad(deg, min, seconds, sign),
                       dso::deg2rad(angle_degrees)));

    // radians to hexicondal degrees and back
    angle_radians = generate_random_double(-dso::D2PI, dso::D2PI);
    dso::rad2hexd(angle_radians, deg, min, seconds, sign);
    assert(approxEqual(dso::hexd2rad(deg, min, seconds, sign), angle_radians));
    assert(approxEqual(dso::hexd2decd(deg, min, seconds, sign),
                       dso::rad2deg(angle_radians)));
  }

  printf("\n");
  return 0;
}
