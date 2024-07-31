#include "transformations.hpp"
#include "units.hpp"
#include <cstdio>

using namespace dso;
const double Re = 6371e3;

int main() {

  return 0;

  Eigen::Matrix<double, 3, 1> maxdiffs = Eigen::Matrix<double, 3, 1>::Zero();
  GeodeticCrd mhgt,mlat,mlon;

  /* first check internal precision with
   * -180 < lon < 180 and -90 < lat < 90
   */
  double hgt = 0e0;
  while (hgt < 10000e0) {
    double lon = -DPI;
    while (lon < DPI) {
      double lat = -DPI / 2e0;
      while (lat < DPI / 2e0) {
        /* geodetic coordinates (φ,λ,h) = (lat, lon, 0) */
        GeodeticCrd s;
        s.lat() = lat;
        s.lon() = lon;
        s.hgt() = hgt;
        /* geodetic to cartesian */
        const auto crt = geodetic2cartesian<ellipsoid::wgs84>(s);
        /* cartesian back to geodetic */
        const auto geo = cartesian2geodetic<ellipsoid::wgs84>(crt);
        /* check diffs */
        if (std::abs(geo.hgt() - s.hgt()) > maxdiffs(0)) {
          maxdiffs(0) = std::abs(geo.hgt() - s.hgt());
          mhgt = s;
        }
        if (std::abs(geo.lat() - s.lat()) > maxdiffs(1)) {
          maxdiffs(1) = std::abs(geo.lat() - s.lat());
          mlat = s;
        }
        if (std::abs(geo.lon() - s.lon()) > maxdiffs(2)) {
          maxdiffs(2) = std::abs(geo.lon() - s.lon());
          mlon = s;
        }
        /* augment latitude */
        lat += 1e-2;
      }
      /* augment lonitude */
      lon += 1e-1;
    }
    /* augment height */
    hgt += 3e0;
  }

  /* report max differences */
  printf("Test : Geodetic -> Cartesian -> Geodetic\n");
  printf("Max height diff = %.2e[m]   at (%.3f[deg], %.3f[deg], %.3f[m])\n",
         maxdiffs(0), rad2deg(mhgt.lat()), rad2deg(mhgt.lon()),
         mhgt.hgt());
  printf("Max lat    diff = %.2e[sec] at (%.3f[deg], %.3f[deg], %.3f[m])\n",
         rad2sec(maxdiffs(1)), rad2deg(mlat.lat()),
         rad2deg(mlat.lon()), mlat.hgt());
  printf("Max lon    diff = %.2e[sec] at (%.3f[deg], %.3f[deg], %.3f[m])\n",
         rad2sec(maxdiffs(2)), rad2deg(mlon.lat()),
         rad2deg(mlon.lon()), mlon.hgt());

  return 0;
}
