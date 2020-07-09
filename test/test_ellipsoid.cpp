#include <iostream>
#include <cassert>
#include <cstring>
#include "ellipsoid.hpp"

using ngpt::ellipsoid;

/*
 * References:
 * [1] H. Moritz, GEODETIC REFERENCE SYSTEM 1980, 
 * https://geodesy.geology.ohio-state.edu/course/refpapers/00740128.pdf
 */

int main()
{

  // here is one way to get ellipsoid (geometric) parameters:
  assert(ngpt::ellipsoid_traits<ellipsoid::grs80>::a == 6378137e0);
  static_assert(ngpt::ellipsoid_traits<ellipsoid::grs80>::a == 6378137e0);
  assert(ngpt::ellipsoid_traits<ellipsoid::wgs84>::f == 1e0/298.257223563e0);
  static_assert(ngpt::ellipsoid_traits<ellipsoid::wgs84>::f == 1e0/298.257223563e0);
  assert(!std::strcmp(ngpt::ellipsoid_traits<ellipsoid::pz90>::n, "PZ90"));

  // for grs80, according to [1] the squared eccentricity should be:
  static_assert(std::abs(ngpt::eccentricity_squared<ellipsoid::grs80>()
                          - .00669438002290)<1e-15);
  // ... and the semi-minor axis is:
  static_assert(std::abs(ngpt::semi_minor<ellipsoid::grs80>() 
                          - 6356752.3141e0)<1e-4);
  // linear eccentricity 
#if defined(__GNUC__) && !defined(__llvm__)
  static_assert(std::abs(ngpt::linear_eccentricity<ellipsoid::grs80>() 
                          - 521854.0097e0)<1e-4);
#else
  assert(std::abs(ngpt::linear_eccentricity<ellipsoid::grs80>() 
                          - 521854.0097e0)<1e-4);
#endif

  return 0;
}
