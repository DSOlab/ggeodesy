#include "geodesy.hpp"
#include "units.hpp"
#include "test_help.hpp"
#include <cassert>
#include <iostream>
#include <random>

using ngpt::normalize_angle;

int main() {
  
  double rand_d, theta;
  double max_sin_dif = std::numeric_limits<double>::min();
  double max_cos_dif = std::numeric_limits<double>::min();
  double mean_sin_dif = 0e0, mean_cos_dif = 0e0;
  std::size_t sin_it=0, cos_it=0;

  /* Normalize in range [0, 2pi] */
  /* ------------------------------------------------------------------------*/
  for (int i=0; i<100000; ++i) {
    rand_d = generate_random_double(-10e0*ngpt::D2PI, 10e0*ngpt::D2PI);
    theta = normalize_angle(rand_d);
    if (!approxEqual(std::sin(theta), std::sin(rand_d))) {
      double dif = std::abs(std::sin(theta)-std::sin(rand_d));
      if (dif > max_sin_dif) max_sin_dif = dif;
      mean_sin_dif = (dif + sin_it*mean_sin_dif) / static_cast<double>(sin_it+1);
      ++sin_it;
    }
    if (!approxEqual(std::cos(theta), std::cos(rand_d))) {
      double dif = std::abs(std::cos(theta)-std::cos(rand_d));
      if (dif > max_cos_dif) max_cos_dif = dif;
      mean_cos_dif = (dif + cos_it*mean_cos_dif) / static_cast<double>(cos_it+1);
      ++cos_it;
    }

    // assert(approxEqual(std::sin(theta), std::sin(rand_d)));
    // assert(approxEqual(std::cos(theta), std::cos(rand_d)));
  }
  printf("Normalizing in range [0, 2*pi]\n");
  printf("Max  cos difference: %+20.17f\n", max_cos_dif);
  printf("Max  sin difference: %+20.17f\n", max_sin_dif);
  printf("Mean cos difference: %+20.17f\n", mean_cos_dif);
  printf("Mean sin difference: %+20.17f\n", mean_sin_dif);
  
  /* test at limits: -2*π */
  rand_d = -ngpt::D2PI;
  theta = normalize_angle(rand_d);
  printf("Special case, normalizing -2pi, diff=%+20.17f\n", std::abs(std::sin(theta)-std::sin(rand_d)));
  printf("                               limit=%+20.17f\n", std::numeric_limits<double>::epsilon() * std::max(std::abs(std::sin(theta)), std::abs(std::sin(rand_d))));
  printf("Special case, normalizing -2pi, diff=%+20.17f\n", std::abs(std::cos(theta)-std::cos(rand_d)));
  printf("                               limit=%+20.17f\n", std::numeric_limits<double>::epsilon() * std::max(std::abs(std::cos(theta)), std::abs(std::cos(rand_d))));
  //assert(approxEqual(std::sin(theta), std::sin(rand_d)));
  //assert(approxEqual(std::cos(theta), std::cos(rand_d)));
  /* test at limits: 2*π */
  rand_d = ngpt::D2PI;
  theta = normalize_angle(rand_d);
  printf("Special case, normalizing 2pi,  diff=%+20.17f\n",  std::abs(std::sin(theta)-std::sin(rand_d)));
  printf("                                limit=%+20.17f\n", std::numeric_limits<double>::epsilon() * std::max(std::abs(std::sin(theta)), std::abs(std::sin(rand_d))));
  printf("Special case, normalizing 2pi,  diff=%+20.17f\n",  std::abs(std::cos(theta)-std::cos(rand_d)));
  printf("                                limit=%+20.17f\n", std::numeric_limits<double>::epsilon() * std::max(std::abs(std::cos(theta)), std::abs(std::cos(rand_d))));
  //assert(approxEqual(std::sin(theta), std::sin(rand_d)));
  //assert(approxEqual(std::cos(theta), std::cos(rand_d)));
  /* test at limits: 0e0 */
  rand_d = 0e0;
  theta = normalize_angle(rand_d);
  printf("Special case, normalizing 0e0,  diff=%+20.17f\n", std::abs(std::sin(theta)-std::sin(rand_d)));
  printf("                               limit=%+20.17f\n", std::numeric_limits<double>::epsilon() * std::max(std::abs(std::sin(theta)), std::abs(std::sin(rand_d))));
  printf("Special case, normalizing 0e0,  diff=%+20.17f\n", std::abs(std::cos(theta)-std::cos(rand_d)));
  printf("                               limit=%+20.17f\n", std::numeric_limits<double>::epsilon() * std::max(std::abs(std::cos(theta)), std::abs(std::cos(rand_d))));
  //assert(approxEqual(std::sin(theta), std::sin(rand_d)));
  //assert(approxEqual(std::cos(theta), std::cos(rand_d)));

  /* Normalize in range [-pi, pi] */
  /* ------------------------------------------------------------------------
  for (int i=0; i<100000; ++i) {
    rand_d = generate_random_double(-10e0*ngpt::D2PI, 10e0*ngpt::D2PI);
    theta = normalize_angle(rand_d, -ngpt::DPI, ngpt::DPI);
    assert(std::abs(std::sin(theta)-std::sin(rand_d))<ErrLimit);
    assert(std::abs(std::cos(theta)-std::cos(rand_d))<ErrLimit);
  }
  / test at limits: -2*π /
  rand_d = -ngpt::D2PI;
  theta = normalize_angle(rand_d, -ngpt::DPI, ngpt::DPI);
  assert(std::abs(std::sin(theta)-std::sin(rand_d))<ErrLimit);
  assert(std::abs(std::cos(theta)-std::cos(rand_d))<ErrLimit);
  / test at limits: -π /
  rand_d = -ngpt::DPI;
  theta = normalize_angle(rand_d, -ngpt::DPI, ngpt::DPI);
  assert(std::abs(std::sin(theta)-std::sin(rand_d))<ErrLimit);
  assert(std::abs(std::cos(theta)-std::cos(rand_d))<ErrLimit);
  / test at limits: π /
  rand_d = ngpt::DPI;
  theta = normalize_angle(rand_d, -ngpt::DPI, ngpt::DPI);
  assert(std::abs(std::sin(theta)-std::sin(rand_d))<ErrLimit);
  assert(std::abs(std::cos(theta)-std::cos(rand_d))<ErrLimit);
  / test at limits: 2*π /
  rand_d = ngpt::D2PI;
  theta = normalize_angle(rand_d, -ngpt::DPI, ngpt::DPI);
  assert(std::abs(std::sin(theta)-std::sin(rand_d))<ErrLimit);
  assert(std::abs(std::cos(theta)-std::cos(rand_d))<ErrLimit);
  / test at limits: 0e0 /
  rand_d = 0e0;
  theta = normalize_angle(rand_d, -ngpt::DPI, ngpt::DPI);
  assert(std::abs(std::sin(theta)-std::sin(rand_d))<ErrLimit);
  assert(std::abs(std::cos(theta)-std::cos(rand_d))<ErrLimit);*/

  return 0;
}
