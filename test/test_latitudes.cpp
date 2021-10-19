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

  /*
  double angle_radians, lat;
  const double a = dso::ellipsoid_traits<dso::ellipsoid::grs80>::a;
  const double f = dso::ellipsoid_traits<dso::ellipsoid::grs80>::f;
  */

  // check that the geocentric and geodetic latitudes are the same at the
  // equator and the poles.
  assert(
      approxEqual(dso::geocentric_latitude<dso::ellipsoid::grs80>(0e0), 0e0));
  assert(approxEqual(
      dso::geocentric_latitude<dso::ellipsoid::grs80>(dso::DPI / 2e0),
      dso::DPI / 2e0));
  dso::Ellipsoid Ell(dso::ellipsoid::grs80);
  assert(
      approxEqual(Ell.geocentric_latitude(-dso::DPI / 2e0), -dso::DPI / 2e0));

  // check that the reduced and geodetic latitudes are the same at the
  // equator and the poles.
  assert(approxEqual(dso::reduced_latitude<dso::ellipsoid::grs80>(0e0), 0e0));
  assert(approxEqual(
      dso::reduced_latitude<dso::ellipsoid::grs80>(dso::DPI / 2e0),
      dso::DPI / 2e0));
  assert(approxEqual(Ell.reduced_latitude(-dso::DPI / 2e0), -dso::DPI / 2e0));

  return 0;
}
