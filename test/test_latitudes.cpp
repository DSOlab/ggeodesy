#include "ellipsoid.hpp"
#include "geoconst.hpp"
#include "test_help.hpp"
#include <cassert>
#include <iostream>

/*
 * The geodetic and geocentric latitudes are equal at the equator and at the
 * poles but at other latitudes they differ by a few minutes of arc. Taking
 * the value of the squared eccentricity as 0.0067 (it depends on the choice
 * of ellipsoid) the maximum difference of φ-θ  may be shown to be about 11.5
 * minutes of arc at a geodetic latitude of approximately 45° 6′
 * see https://en.wikipedia.org/wiki/Latitude
 */

int main() {

  double angle_radians, lat;
  const double a = ngpt::ellipsoid_traits<ngpt::ellipsoid::grs80>::a;
  const double f = ngpt::ellipsoid_traits<ngpt::ellipsoid::grs80>::f;

  // check that the geocentric and geodetic latitudes are the same at the
  // equator and the poles.
  assert(
      approxEqual(ngpt::geocentric_latitude<ngpt::ellipsoid::grs80>(0e0), 0e0));
  assert(approxEqual(
      ngpt::geocentric_latitude<ngpt::ellipsoid::grs80>(ngpt::DPI / 2e0),
      ngpt::DPI / 2e0));
  ngpt::Ellipsoid Ell(ngpt::ellipsoid::grs80);
  assert(
      approxEqual(Ell.geocentric_latitude(-ngpt::DPI / 2e0), -ngpt::DPI / 2e0));

  // check that the reduced and geodetic latitudes are the same at the
  // equator and the poles.
  assert(approxEqual(ngpt::reduced_latitude<ngpt::ellipsoid::grs80>(0e0), 0e0));
  assert(approxEqual(
      ngpt::reduced_latitude<ngpt::ellipsoid::grs80>(ngpt::DPI / 2e0),
      ngpt::DPI / 2e0));
  assert(approxEqual(Ell.reduced_latitude(-ngpt::DPI / 2e0), -ngpt::DPI / 2e0));

  return 0;
}
