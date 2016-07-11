#include "ellipsoid.hpp"
#include "car2ell.hpp"
#include "ell2car.hpp"
#include "car2top.hpp"
#include "vincenty.hpp"
#include "geodesy.hpp"

#include <stdio.h>
#include <cmath>
#include <cassert>

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

    // --------------------------------------------------------------------- //
    //                     TEST INVERSE VINCENTY ALGORITHM
    // --------------------------------------------------------------------- //

    double S, a_for, a_bac, a_sec;
    int a_deg, a_min;
    // try the (inverse) vincenty algorithm; test data/output is taken from
    // http://www.ga.gov.au/geodesy/datums/vincenty_inverse.jsp
    p1.x = -(37.0 + 57/60.0e0 + 3.72030/3600.0e0);
    p1.y = 144.0  + 25/60.0e0 + 29.52440/3600.0e0;
    p2.x = -(37.0 + 39/60.0e0 + 10.15610/3600.0e0);
    p2.y = 143.0  + 55/60.0e0 + 35.38390/3600.0e0;
    p1.x *= DPI/180.0;
    p1.y *= DPI/180.0;
    p2.x *= DPI/180.0;
    p2.y *= DPI/180.0;
    // Results should be:
    // forward azimouth:  306° 52' 5.37
    // backward azimouth: 127° 10' 25.07
    // distance         : 54972.271 m
    S = inverse_vincenty(p1.x, p1.y, p2.x, p2.y, a_for, a_bac, 1e-12);
    
    assert( std::abs(S-54972.271) < 1e-3 ); // millimeter accuracy
    assert( std::abs(a_for-hexd2rad(306, 52, 5.37)) < 3*1e-5 ); // .1 sec accuracy
    assert( std::abs(a_bac-hexd2rad(127, 10,25.07)) < 3*1e-5 ); // .1 sec accuracy
    
    printf("Distance of the ellipsoidal: %+15.5f m\n", S);
    rad2hexd(a_for, a_deg, a_min, a_sec);
    printf("Forward Azimouth           : %+3d %2d %8.5f (or %15.10f rad)\n", a_deg, a_min, a_sec, a_for);
    rad2hexd(a_bac, a_deg, a_min, a_sec);
    printf("Backward Azimouth          : %+3d %2d %8.5f (or %15.10f rad)\n", a_deg, a_min, a_sec, a_bac);
    
    // --------------------------------------------------------------------- //
    //                     TEST VINCENTY ALGORITHM
    // --------------------------------------------------------------------- //
    // Do the inverse from above (so results should be the same, going back
    // to p2)
    double new_lon, new_lat, new_az;
    new_az = direct_vincenty(p1.x, p1.y, a_for, S, new_lat, new_lon, 1e-12);
    //assert(std::abs(new_lon-p2.x) < 3*1e-7 ); // .001 sec accuracy
    //assert(std::abs(new_lat-p2.y) < 3*1e-7 ); // .001 sec accuracy
    //assert(std::abs(new_az-a_bac) < 3*1e-7 ); // .001 sec accuracy

    rad2hexd(new_lon, a_deg, a_min, a_sec);
    printf("Longtitude:  %+3d %2d %8.5f\n", a_deg, a_min, a_sec);
    rad2hexd(new_lat, a_deg, a_min, a_sec);
    printf("Latitude  :  %+3d %2d %8.5f\n", a_deg, a_min, a_sec);
    rad2hexd(new_az, a_deg, a_min, a_sec);
    printf("Azimouth  :  %+3d %2d %8.5f\n", a_deg, a_min, a_sec);

    return 0;
}
