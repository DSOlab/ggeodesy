#include "geodesy.hpp"
#include <cassert>
#include <iostream>
#include <random>

using ngpt::normalize_angle;

double generate_random_double(double lower_bound, double upper_bound) noexcept {
  std::random_device rd;
  std::mt19937 e2(rd());
  std::uniform_real_distribution<double> unif(lower_bound, upper_bound);
  return unif(e2);
}

constexpr double ErrLimit = 1e-13;

int main() {
  
  double rand_d, theta;

  /* Normalize in range [0, 2pi] */
  /* ------------------------------------------------------------------------*/
  for (int i=0; i<100000; ++i) {
    rand_d = generate_random_double(-10e0*ngpt::D2PI, 10e0*ngpt::D2PI);
    theta = normalize_angle(rand_d);
    assert(std::abs(std::sin(theta)-std::sin(rand_d))<ErrLimit);
    assert(std::abs(std::cos(theta)-std::cos(rand_d))<ErrLimit);
  }
  /* test at limits: -2*π */
  rand_d = -ngpt::D2PI;
  theta = normalize_angle(rand_d);
  assert(std::abs(std::sin(theta)-std::sin(rand_d))<ErrLimit);
  assert(std::abs(std::cos(theta)-std::cos(rand_d))<ErrLimit);
  /* test at limits: 2*π */
  rand_d = ngpt::D2PI;
  theta = normalize_angle(rand_d);
  assert(std::abs(std::sin(theta)-std::sin(rand_d))<ErrLimit);
  assert(std::abs(std::cos(theta)-std::cos(rand_d))<ErrLimit);
  /* test at limits: 0e0 */
  rand_d = 0e0;
  theta = normalize_angle(rand_d);
  assert(std::abs(std::sin(theta)-std::sin(rand_d))<ErrLimit);
  assert(std::abs(std::cos(theta)-std::cos(rand_d))<ErrLimit);
  
  /* Normalize in range [-pi, pi] */
  /* ------------------------------------------------------------------------*/
  for (int i=0; i<100000; ++i) {
    rand_d = generate_random_double(-10e0*ngpt::D2PI, 10e0*ngpt::D2PI);
    theta = normalize_angle(rand_d, -ngpt::DPI, ngpt::DPI);
    assert(std::abs(std::sin(theta)-std::sin(rand_d))<ErrLimit);
    assert(std::abs(std::cos(theta)-std::cos(rand_d))<ErrLimit);
  }
  /* test at limits: -2*π */
  rand_d = -ngpt::D2PI;
  theta = normalize_angle(rand_d, -ngpt::DPI, ngpt::DPI);
  assert(std::abs(std::sin(theta)-std::sin(rand_d))<ErrLimit);
  assert(std::abs(std::cos(theta)-std::cos(rand_d))<ErrLimit);
  /* test at limits: -π */
  rand_d = -ngpt::DPI;
  theta = normalize_angle(rand_d, -ngpt::DPI, ngpt::DPI);
  assert(std::abs(std::sin(theta)-std::sin(rand_d))<ErrLimit);
  assert(std::abs(std::cos(theta)-std::cos(rand_d))<ErrLimit);
  /* test at limits: π */
  rand_d = ngpt::DPI;
  theta = normalize_angle(rand_d, -ngpt::DPI, ngpt::DPI);
  assert(std::abs(std::sin(theta)-std::sin(rand_d))<ErrLimit);
  assert(std::abs(std::cos(theta)-std::cos(rand_d))<ErrLimit);
  /* test at limits: 2*π */
  rand_d = ngpt::D2PI;
  theta = normalize_angle(rand_d, -ngpt::DPI, ngpt::DPI);
  assert(std::abs(std::sin(theta)-std::sin(rand_d))<ErrLimit);
  assert(std::abs(std::cos(theta)-std::cos(rand_d))<ErrLimit);
  /* test at limits: 0e0 */
  rand_d = 0e0;
  theta = normalize_angle(rand_d, -ngpt::DPI, ngpt::DPI);
  assert(std::abs(std::sin(theta)-std::sin(rand_d))<ErrLimit);
  assert(std::abs(std::cos(theta)-std::cos(rand_d))<ErrLimit);

  return 0;
}
