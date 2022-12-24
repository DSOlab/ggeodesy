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
 * Start with an angle in radians, transform it to degrees and then back to
 * radians. See the difference
 */
int main() {
  std::uniform_real_distribution<double> unif(0e0,1e-5);
  std::default_random_engine re;

  double maxdif = std::numeric_limits<double>::min();
  double maxdifangle=0e0;
  double d;
  int diffound = false;

  for (double a=-8*M_PI; a<8*M_PI; a += 1e-7) {
    const double angle = a + unif(re);
    // to degrees ...
    const double deg = dso::rad2deg<double>(angle);
    // and back to radians ...
    const double rad = dso::deg2rad<double>(deg);
    // are they the same number ?
    if (!ApproxEqual(angle, rad)) {
      ++diffound;
      if (maxdif < (d = std::abs(rad - angle))) {
        maxdif = d;
        maxdifangle = angle;
      }
    }
  }

  if (diffound) {
    printf("Max difference is %.9e at angle %.15e[rad] or %.15e[deg]\n", maxdif,
           maxdifangle, dso::rad2deg(maxdifangle));
    fprintf(stderr, "Transformation rad->deg->rad NOT within machine "
                    "precision. I will be failing now ...\n");
    assert(false);
  } else {
    printf("All transformations rad->deg->rad produce results within machine "
           "precision\n");
  }

  return 0;
}
