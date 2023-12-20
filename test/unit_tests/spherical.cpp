#include "geodesy.hpp"
#include "units.hpp"
#include <cstdio>

const double Re = 6371e3;

int main() {

  Eigen::Matrix<double, 3, 1> maxdiffs = Eigen::Matrix<double, 3, 1>::Zero();

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
      if (std::abs(sph.r() - s.r()) > maxdiffs(0))
        maxdiffs(0) = std::abs(sph.r() - s.r());
      if (std::abs(sph.lat() - s.lat()) > maxdiffs(1))
        maxdiffs(1) = std::abs(sph.lat() - s.lat());
      if (std::abs(sph.lon() - s.lon()) > maxdiffs(2))
        maxdiffs(2) = std::abs(sph.lon() - s.lon());
      /* augment latitude */
      lat += 1e-3;
    }
    /* augment lonitude */
    lon += 1e-3;
  }

  /* report max differences */
  printf("Max diffs: %.3e[sec] %.3e[sec] %.3e[mm]\n", dso::rad2sec(maxdiffs(1)),
         dso::rad2sec(maxdiffs(2)), maxdiffs(0) * 1e3);

  /* second, report results for ranges
   * -360 < lon < 360 and -360 < lat < 360
   */
  lon = -2 * dso::DPI;
  while (lon < 2 * dso::DPI) {
    double lat = -2 * dso::DPI;
    while (lat < 2 * dso::DPI) {
      /* spherical coordinates (r,φ,λ) = (Re, lat, lon) */
      dso::SphericalCrd s;
      s.r() = Re;
      s.lat() = lat;
      s.lon() = lon;
      /* spherical to cartesian */
      const auto crt = dso::spherical2cartesian(s);
      printf("lat=%+.20e lon=%+.20e R=%.20e X=%+.17e Y=%.17e Z=%.17e\n",
             s.lat(), s.lon(), s.r(), crt.x(), crt.y(), crt.z());
      /* augment latitude */
      lat += 1e-2;
    }
    /* augment lonitude */
    lon += 1e-2;
  }

  return 0;
}
