#include "ellipsoid.hpp"
#include <cassert>
#include <cstring>
#include <iostream>

using ngpt::ellipsoid;
using ngpt::Ellipsoid;

/*
 * References:
 * [1] H. Moritz, GEODETIC REFERENCE SYSTEM 1980,
 * https://geodesy.geology.ohio-state.edu/course/refpapers/00740128.pdf
 */

int main() {

  // Using the enum (class) ellipsoid
  // here is one way to get ellipsoid (geometric) parameters:
  assert(ngpt::ellipsoid_traits<ellipsoid::grs80>::a == 6378137e0);
  static_assert(ngpt::ellipsoid_traits<ellipsoid::grs80>::a == 6378137e0);
  assert(ngpt::ellipsoid_traits<ellipsoid::wgs84>::f == 1e0 / 298.257223563e0);
  static_assert(ngpt::ellipsoid_traits<ellipsoid::wgs84>::f ==
                1e0 / 298.257223563e0);
  assert(!std::strcmp(ngpt::ellipsoid_traits<ellipsoid::pz90>::n, "PZ90"));

  // for grs80, according to [1] the squared eccentricity should be:
  static_assert(std::abs(ngpt::eccentricity_squared<ellipsoid::grs80>() -
                         .00669438002290) < 1e-15);
  // ... and the semi-minor axis is:
  static_assert(
      std::abs(ngpt::semi_minor<ellipsoid::grs80>() - 6356752.3141e0) < 1e-4);
  // linear eccentricity
#if defined(__GNUC__) && !defined(__llvm__)
  static_assert(std::abs(ngpt::linear_eccentricity<ellipsoid::grs80>() -
                         521854.0097e0) < 1e-4);
#else
  assert(std::abs(ngpt::linear_eccentricity<ellipsoid::grs80>() -
                  521854.0097e0) < 1e-4);
#endif
  // ... and the polar radius of curvature is:
  static_assert(std::abs(ngpt::polar_radius_of_curvature<ellipsoid::grs80>() -
                         6399593.6259) < 1e-4);

  // using the Ellipsoid class
  // declare Ellipsoid (class) instances
  constexpr auto wgs84 = Ellipsoid(ellipsoid::wgs84);
  constexpr auto grs80 = Ellipsoid(ngpt::ellipsoid_traits<ellipsoid::grs80>::a,
                                   ngpt::ellipsoid_traits<ellipsoid::grs80>::f);
  constexpr auto pz90 = Ellipsoid(ellipsoid::pz90);

  static_assert(wgs84.eccentricity_squared() ==
                ngpt::eccentricity_squared<ellipsoid::wgs84>());
  static_assert(grs80.semi_minor() == ngpt::semi_minor<ellipsoid::grs80>());
  static_assert(std::abs(pz90.eccentricity_squared() - 0.0066943662) < 1e-9);

  return 0;
}
