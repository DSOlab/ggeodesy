#include "car2ell.hpp"
#include "car2top.hpp"
#include "geodesy.hpp"
#include "units.hpp"
#include "ellipsoid.hpp"
#include <cassert>
#include <stdio.h>

struct Point {
  double x, y, z;
};

// distance between two points
double distance(const Point &p1, const Point &p2) noexcept {
  return std::sqrt((p1.x - p2.x) * (p1.x - p2.x) +
                   (p1.y - p2.y) * (p1.y - p2.y) +
                   (p1.z - p2.z) * (p1.z - p2.z));
}

// norm of vector
double distance(const Point &p1) noexcept {
  return std::sqrt(p1.x * p1.x + p1.y * p1.y + p1.z * p1.z);
}

// precision for comparing two doubles for equality
constexpr double prec{1e-10};
constexpr double min_prec{1e-4};
constexpr double max_prec{1e-15};

using namespace ngpt;

int main() {
  /// Points on ellipsoid.
  Point p1{4595220.02261, 2039434.13622, 3912625.96044};
  Point p2{4595220.02261 + 5.12345, 2039434.13622 + 5.12345,
           3912625.96044 + 5.12345};
  // Used to check precision.
  double p;

  // --------------------------------------------------------------------- //
  //                     TEST COORDINATE TRANSFORMATIONS
  // --------------------------------------------------------------------- //
  printf("\n> TEST COORDINATE TRANSFORMATIONS\n");
  printf("------------------------------------\n");
  printf("\n\nWe start with the following two points (cartesian):");
  printf("\n| x |   |%+15.5f| | x |   |%+15.5f|", p1.x, p2.x);
  printf("\n| y | = |%+15.5f| | y | = |%+15.5f|", p1.y, p2.y);
  printf("\n| z |A  |%+15.5f| | z |B  |%+15.5f|", p1.z, p2.z);

  printf("\n\nCartesian and Topocentric vectors, between Pa and Pb (i.e Pb-Pa) "
         "are:");
  //  P3 is the [dx,dy,dz] vector
  Point p3{p2.x - p1.x, p2.y - p1.y, p2.z - p1.z};
  //  P4 will hold the topocentric vector, [north, east, up]
  Point p4;
  //  Transform P2-P1 to a topocentric vector around P1; store results in P4
  //  Note that here we pass both points (p1 and P2) to the function.
  car2top<ellipsoid::grs80>(p1.x, p1.y, p1.z, p2.x, p2.y, p2.z, p4.x, p4.y,
                            p4.z);
  //  However, since we already have the [dx,dy,dz] vector, we could transform
  //+ it directly to topocentric; i.e. instead of using:
  //+ car2top<ellipsoid::grs80>(p1.x, p1.y, p1.z, p2.x, p2.y, p2.z, p4.x, p4.y,
  // p4.z);
  //+ we could have used:
  Point p4_2;
  dcar2top<ellipsoid::grs80>(p1.x, p1.y, p1.z, p3.x, p3.y, p3.z, p4_2.x, p4_2.y,
                             p4_2.z);
  // Now verify that both functions give the same results
  assert(std::abs(p4_2.x - p4.x) < prec && std::abs(p4_2.y - p4.y) < prec &&
         std::abs(p4_2.z - p4.z) < prec);
  printf("\n| δx |   |%+15.10f|  | n |   |%+15.10f|", p3.x, p4.x);
  printf("\n| δy | = |%+15.10f|  | e | = |%+15.10f|", p3.y, p4.y);
  printf("\n| δz |AB |%+15.10f|  | u |AB |%+15.10f|", p3.z, p4.z);
  //  Check to what precision we get the same results
  for (p = max_prec; p >= min_prec; p *= 1e-1) {
    if (!(std::abs(p4_2.x - p4.x) < prec && std::abs(p4_2.y - p4.y) < prec &&
          std::abs(p4_2.z - p4.z) < prec))
      break;
  }
  printf("\nInternal precision of the two car2top implementations, is "
         "approximately %e meters",
         p);

  //  The distance between two points, is the same in both geocentric and
  //+ local reference systems.
  double d1, d2, d3;
  //  Distance between P1 and P2, using geocentric cartesian coordinates
  d1 = distance(p1, p2);
  //  Norm of p3 vector (i.e. the [dx,dy,dz] vector)
  d2 = distance(p3);
  //  Noerm of p4 vector (i.e. the [north, east, up] vector
  d3 = distance(p4);
  printf("\n\nDistance between points p1 and p2, calculated from cartesian");
  printf("\ncomponents, is                        : %15.7f", d1);
  printf("\nNorm of the cartesian diff. vector is : %15.7f", d2);
  printf("\nNorm of the topocentric vector is     : %15.7f", d3);
  //  Verify the results
  assert(std::abs(d1 - d2) < prec && std::abs(d2 - d3) < prec);

  //  From the local, topocentric vector p3 = [north, east, up], we can
  //  calculate
  //+ the distance, azimouth and zenith angle between the two points. Obviously
  //+ the distance should match the ones we computed above (d1, d2 and d3).
  printf("\n\nFrom the topocentric [N, E, U] vector, we can calculate the");
  printf("\ndistance, azimouth and zenith angle between the two points:");
  double dist, azi, zen;
  top2daz(p4.x, p4.y, p4.z, dist, azi, zen);
  printf("\ndistance   = %15.7f (m)", dist);
  printf("\nazimouth   = %15.10f deg.", rad2deg(azi));
  printf("\nzenith ang.= %15.10f deg.", rad2deg(zen));
  //  Verify results
  assert(std::abs(dist - d1) < prec);

  //  We can also transform a local topocentric vector to a geocentric cartesian
  //+ one. For this, we will need the function top2car, and the ellipsoidal
  //+ coordinates of the reference point (i.e. P1).
  printf("\n\nFrom the topocentric vector we can go back to cartesian, if we");
  printf("\nknow the coordinates (actually latitude and longtitude of the");
  printf("\nreference point. Hence, actually compute the cartesian components");
  printf("\nof the second, ending point");
  Point p5, // ellipsoidal coordinates of p1
      p6;   // cartesian coordinates of p2
  //  Get ellipsoidal coordinates (lat, lon, height) of P1 and store them to
  //+ P5
  car2ell<ellipsoid::grs80>(p1.x, p1.y, p1.z, p5.x, p5.y, p5.z);
  //  Transform the topocentric vector p3 = [north, east, up] to a geocentric
  //+ cartesian one, with P1 as reference point fot the local coordinate sys.
  top2car(p4.x, p4.y, p4.z, p5.x, p5.y, p6.x, p6.y, p6.z);
  printf("\n| φ |   |%+15.10f| | δx |    |%+15.10f|", rad2deg(p5.x), p6.x);
  printf("\n| λ | = |%+15.10f| | δy | =  |%+15.10f|", rad2deg(p5.y), p6.y);
  printf("\n| h |A  |%+14.5f | | δz |AB  |%+15.10f|", p5.z, p6.z);
  //  Verify that the results are the same as the original Dx vector, i.e. P3
  assert(std::abs(p6.x - p3.x) < prec && std::abs(p6.y - p3.y) < prec &&
         std::abs(p6.z - p3.z) < prec);
  //  Check to what precision we get the same results
  for (p = max_prec; p >= min_prec; p *= 1e-1) {
    if (!(std::abs(p6.x - p6.x) < prec && std::abs(p6.y - p3.y) < prec &&
          std::abs(p6.z - p6.z) < prec))
      break;
  }
  printf("\nInternal precision of top2car, is approximately %e meters", p);

  printf("\n");
  return 0;
}
