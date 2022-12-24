#include "units.hpp"
#include <cassert>
#include <cstdio>
#include <limits>
#include <cmath>
#include <random>

bool ApproxEqual(double a, double b,
                 double tolerance = std::numeric_limits<double>::epsilon()) {
  double diff = std::abs(a - b);
  if (diff <= tolerance)
    return true;

  if (diff < std::max(std::abs(a), std::abs(b)) * tolerance)
    return true;

  return false;
}

/*
 * Start with an angle in degrees, transform it to radians and then back to
 * degrees. See the difference
 */
int main() {
  std::uniform_real_distribution<double> unif(0e0,1e-4);
  std::default_random_engine re;

  double maxdif = std::numeric_limits<double>::min();
  double maxdifangle=0e0;
  double d;
  int diffound = false;

  for (double a=-4*360e0; a<4*360e0; a += 1e-4) {
    const double angle = a + unif(re);
    // to radians ...
    const double rad = dso::deg2rad<double>(angle);
    // and back to degrees ...
    const double deg = dso::rad2deg<double>(rad);
    // are they the same number ?
    if (!ApproxEqual(angle, deg)) {
      ++diffound;
      if (maxdif < (d = std::abs(deg - angle))) {
        maxdif = d;
        maxdifangle = angle;
      }
    }
  }

  if (diffound) {
    printf("Max difference is %.9e at angle %.15e[deg] or %.15e[rad]\n", maxdif,
           maxdifangle, dso::deg2rad(maxdifangle));
    fprintf(stderr, "Transformation deg->rad->deg NOT within machine "
                    "precision. I will be failing now ...\n");
    assert(false);
  } else {
    printf("All transformations deg->rad->deg produce results within machine "
           "precision\n");
  }

  return 0;
}
