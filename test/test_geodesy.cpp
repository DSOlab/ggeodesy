#include "ellipsoid.hpp"
#include "car2ell.hpp"
#include "ell2car.hpp"
#include "car2top.hpp"
#include "vincenty.hpp"

#include <stdio.h>
#include <cmath>

struct Point { double x,y,z; };

using namespace ngpt;

int main ()
{
    /// a (reference) point on ellipsoid.
    Point p1 {4595220.02261, 2039434.13622, 3912625.96044};
    Point p2, p3;

    /// cartesian to ellipsoidal; p2 holds the ellipsoidal coordinates of p1.
    car2ell<ellipsoid::grs80>(p1.x,p1.y,p1.z,p2.x,p2.y,p2.z);
    /// ellipsoidal to cartesian; p3 hold the cartesian coordinates of p2 (i.e.
    /// p3 should be the same as p1).
    ell2car<ellipsoid::grs80>(p2.x,p2.y,p2.z,p3.x,p3.y,p3.z);

    /// print the difference between them components of p1 and p3
    printf ("To ellipsoidal and back (i.e. everything should be zero!):");
    printf ("\ndx=%8.5f dy=%8.5F dz=%8.5f\n",
            std::abs(p1.x-p3.x),std::abs(p1.y-p3.y),std::abs(p1.z-p3.z));

    /// compute the topocentric vector between p1 and p3 ...
    car2top<ellipsoid::grs80>(p1.x,p1.y,p1.z,p3.x,p3.y,p3.z,p2.x,p2.y,p2.z);
    /// ... and print the components (should be zero)
    printf ("Topocentric vector (should be zero!)");
    printf ("\ndn=%8.5f de=%8.5f du=%8.5f\n",p2.x,p2.y,p2.z);

    // try the vincenty algorithm; test data/output is taken from
    // http://www.ga.gov.au/geodesy/datums/vincenty_inverse.jsp
    p1.x = -(37.0 + 57/60.0e0 + 3.72030/3600.0e0);
    p1.y = 144.0  + 25/60.0e0 + 29.52440/3600.0e0;
    p2.x = -(37.0 + 39/60.0e0 + 10.15610/3600.0e0);
    p2.y = 143.0  + 55/60.0e0 + 35.38390/3600.0e0;
    p1.x *= DPI/180.0;
    p1.y *= DPI/180.0;
    p2.x *= DPI/180.0;
    p2.y *= DPI/180.0;

    double S, a_for, a_bac;
    S = inverse_vincenty(p1.x, p1.y, p2.x, p2.y, a_for, a_bac, 1e-12);


    return 0;
}
