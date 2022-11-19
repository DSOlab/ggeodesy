#include "geodesy.hpp"
#include "units.hpp"
#include <cstdio>

int main() {
  // angles from [0 to 90]
  double angle = 0e0;
  
  while (angle <= 90e0) {
    const double angrad = dso::deg2rad(angle);
    
    const double geocentric =
        dso::geocentric_latitude<dso::ellipsoid::grs80>(angrad);

    const double reduced = dso::reduced_latitude<dso::ellipsoid::grs80>(angrad);

    const double ellgeo = dso::rad2sec(angrad-geocentric);
    const double redgeo = dso::rad2sec(angrad-reduced);

    printf("%.3f %.12f %.12f %.6f %.6f\n", angle, dso::rad2deg(geocentric),
           dso::rad2deg(reduced), ellgeo, redgeo);

    angle += .05e0;
  }

  return 0;
}
