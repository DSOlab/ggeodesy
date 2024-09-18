#include "transformations.hpp"
#include "units.hpp"
#include <cassert>
#include <cstdio>

using namespace dso;
const double Re = 6371e3;
constexpr const double MAX_DIFF_LAT_ASEC = 1e-8;
constexpr const double MAX_DIFF_LON_ASEC = 5e-11;
constexpr const double MAX_DIFF_HGT_MTRS = 5e-9;

constexpr const double MAX_DIFF_LAT_RAD = sec2rad(MAX_DIFF_LAT_ASEC);
constexpr const double MAX_DIFF_LON_RAD = sec2rad(MAX_DIFF_LON_ASEC);

int main() {

  Eigen::Matrix<double, 3, 1> maxdiffs = Eigen::Matrix<double, 3, 1>::Zero();
  dso::SphericalCrd mlat, mlon, mhgt;

  /* first check internal precision with
   * -180 < lon < 180 and -90 < lat < 90
   */
  double lon = -dso::DPI;
  while (lon < dso::DPI) {
    double lat = -dso::DPI / 2e0;
    while (lat < dso::DPI / 2e0) {
      /* spherical coordinates (r,φ,λ) = (Re, lat, lon) */
      dso::SphericalCrd s;
      s.r() = Re;
      s.lat() = lat;
      s.lon() = lon;
      /* spherical to cartesian */
      const auto crt = dso::spherical2cartesian(s);
      /* cartesian back to spherical */
      const auto sph = dso::cartesian2spherical(crt);
      /* check diffs */
      if (std::abs(sph.r() - s.r()) > maxdiffs(0)) {
        maxdiffs(0) = std::abs(sph.r() - s.r());
#ifndef REPORT_TEST_DIFFS
        assert(maxdiffs(0) < MAX_DIFF_HGT_MTRS);
#endif
        mhgt = s;
      }
      if (std::abs(sph.lat() - s.lat()) > maxdiffs(1)) {
        maxdiffs(1) = std::abs(sph.lat() - s.lat());
#ifndef REPORT_TEST_DIFFS
        assert(maxdiffs(1) < MAX_DIFF_LAT_RAD);
#endif
        mlat = s;
      }
      if (std::abs(sph.lon() - s.lon()) > maxdiffs(2)) {
        maxdiffs(2) = std::abs(sph.lon() - s.lon());
#ifndef REPORT_TEST_DIFFS
        assert(maxdiffs(2) < MAX_DIFF_LON_RAD);
#endif
        mlon = s;
      }
      /* augment latitude */
      lat += 1e-2;
    }
    /* augment lonitude */
    lon += 1e-2;
  }

#ifdef REPORT_TEST_DIFFS
  /* report max differences */
  printf("Test : Spherical -> Cartesian -> Spherical\n");
  printf("Max heihgt diff = %.2e[m]   at (%.3f[deg], %.3f[deg], %.3f[m])\n",
         maxdiffs(0), dso::rad2deg(mhgt.lat()), dso::rad2deg(mhgt.lon()),
         mhgt.r());
  printf("Max lat    diff = %.2e[sec] at (%.3f[deg], %.3f[deg], %.3f[m])\n",
         dso::rad2sec(maxdiffs(1)), dso::rad2deg(mlat.lat()),
         dso::rad2deg(mlat.lon()), mlat.r());
  printf("Max lon    diff = %.2e[sec] at (%.3f[deg], %.3f[deg], %.3f[m])\n",
         dso::rad2sec(maxdiffs(2)), dso::rad2deg(mlon.lat()),
         dso::rad2deg(mlon.lon()), mlon.r());
#endif

  return 0;
}
