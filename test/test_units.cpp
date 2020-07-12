#include "geodesy.hpp"
#include <cassert>
#include <iostream>
#include <random>

double generate_random_double(double lower_bound, double upper_bound) noexcept {
  std::random_device rd;
  std::mt19937 e2(rd());
  std::uniform_real_distribution<double> unif(lower_bound, upper_bound);
  return unif(e2);
}

constexpr double ErrLimit = 1e-13;

int main() {
  double angle_radians, angle_degrees, seconds, theta;
  int deg, min;

  for (int i = 0; i < 150; ++i) {
    // radians to dec. degrees and back
    angle_radians = generate_random_double(-ngpt::D2PI, ngpt::D2PI);
    angle_degrees = ngpt::rad2deg(angle_radians);
    assert(std::abs(ngpt::deg2rad(angle_degrees) - angle_radians) < ErrLimit);

    // dec. degrees to radians and back
    angle_degrees = generate_random_double(-360e0, 360e0);
    angle_radians = ngpt::deg2rad(angle_degrees);
    assert(std::abs(ngpt::rad2deg(angle_radians) - angle_degrees) < ErrLimit);

    // dec. degrees to hexicondal degrees and back
    angle_degrees = generate_random_double(-360e0, 360e0);
    ngpt::decd2hexd(angle_degrees, deg, min, seconds);
    assert(std::abs(ngpt::hexd2decd(deg, min, seconds) - angle_degrees) <
           ErrLimit);
    assert(std::abs(ngpt::hexd2rad(deg, min, seconds) -
                    ngpt::deg2rad(angle_degrees)) < ErrLimit);

    // radians to hexicondal degrees and back
    angle_radians = generate_random_double(-ngpt::D2PI, ngpt::D2PI);
    ngpt::rad2hexd(angle_radians, deg, min, seconds);
    assert(std::abs(ngpt::hexd2rad(deg, min, seconds) - angle_radians) <
           ErrLimit);
    assert(std::abs(ngpt::hexd2decd(deg, min, seconds) -
                    ngpt::rad2deg(angle_radians)) < ErrLimit);

    // check the angle normalization
    angle_radians = generate_random_double(-ngpt::D2PI, ngpt::D2PI);
    theta = ngpt::normalize_angle(angle_radians, 0e0, ngpt::D2PI); // [0, 2π)
    assert(std::abs(std::sin(angle_radians) - std::sin(theta)) < ErrLimit);
    assert(std::abs(std::cos(angle_radians) - std::cos(theta)) < ErrLimit);
    theta =
        ngpt::normalize_angle(angle_radians, -ngpt::DPI, ngpt::DPI); // [-π, π)
    assert(std::abs(std::sin(angle_radians) - std::sin(theta)) < ErrLimit);
    assert(std::abs(std::cos(angle_radians) - std::cos(theta)) < ErrLimit);
  }

  return 0;
}
