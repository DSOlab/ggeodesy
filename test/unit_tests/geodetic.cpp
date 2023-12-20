#include "geodesy.hpp"
#include "units.hpp"
#include <cstdio>

const double Re = 6371e3;
using dso::ellipsoid;

int main() {

  Eigen::Matrix<double, 3, 1> maxdiffs = Eigen::Matrix<double, 3, 1>::Zero();

  /* first check internal precision with
   * -180 < lon < 180 and -90 < lat < 90
   */
  double hgt = 0e0;
  while (hgt < 10000e0) {
    double lon = -dso::DPI;
    while (lon < dso::DPI) {
      double lat = -dso::DPI / 2e0;
      while (lat < dso::DPI / 2e0) {
        /* geodetic coordinates (φ,λ,h) = (lat, lon, 0) */
        dso::GeodeticCrd s;
        s.lat() = lat;
        s.lon() = lon;
        s.hgt() = hgt;
        /* geodetic to cartesian */
        const auto crt = dso::geodetic2cartesian<ellipsoid::wgs84>(s);
        /* cartesian back to geodetic */
        const auto geo = dso::cartesian2geodetic<ellipsoid::wgs84>(crt);
        /* check diffs */
        if (std::abs(geo.hgt() - s.hgt()) > maxdiffs(0))
          maxdiffs(0) = std::abs(geo.hgt() - s.hgt());
        if (std::abs(geo.lat() - s.lat()) > maxdiffs(1))
          maxdiffs(1) = std::abs(geo.lat() - s.lat());
        if (std::abs(geo.lon() - s.lon()) > maxdiffs(2))
          maxdiffs(2) = std::abs(geo.lon() - s.lon());
        /* augment latitude */
        lat += 1e-3;
      }
      /* augment lonitude */
      lon += 1e-3;
    }
    /* augment height */
    hgt += 1e2;
  }

  /* report max differences */
  printf("Max diffs: %.3e[sec] %.3e[sec] %.3e[mm]\n", dso::rad2sec(maxdiffs(1)),
         dso::rad2sec(maxdiffs(2)), maxdiffs(0) * 1e3);

  /* second, report results for ranges
   * -360 < lon < 360 and -360 < lat < 360
   */
  // lon = -2 * dso::DPI;
  // while (lon < 2 * dso::DPI) {
  //   double lat = -2 * dso::DPI;
  //   while (lat < 2 * dso::DPI) {
  //     /* spherical coordinates (r,φ,λ) = (Re, lat, lon) */
  //     dso::SphericalCrd s;
  //     s.r() = Re;
  //     s.lat() = lat;
  //     s.lon() = lon;
  //     /* spherical to cartesian */
  //     const auto crt = dso::spherical2cartesian(s);
  //     printf("lat=%+.20e lon=%+.20e R=%.20e X=%+.17e Y=%.17e Z=%.17e\n",
  //            s.lat(), s.lon(), s.r(), crt.x(), crt.y(), crt.z());
  //     /* augment latitude */
  //     lat += 1e-2;
  //   }
  //   /* augment lonitude */
  //   lon += 1e-2;
  // }

  return 0;
}
