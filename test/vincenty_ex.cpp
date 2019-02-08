#include "geoconst.hpp"
#include "vincenty.hpp"
#include "geodesy.hpp"
#include "trnsfdtls.hpp"

#include <stdio.h>
#include <cmath>
#include <cassert>

struct Point { double x,y,z; };

using namespace ngpt;

int main ()
{
    Point  p1, p2;
    double S,a_for(0e0),a_bac(0e0),a_sec(0e0),a1_sec,a2_sec,b1_sec,b2_sec,
           new_lon,new_lat,new_az;
    int    a_deg,a_min,a1_deg,a1_min,a2_deg,a2_min,b1_deg,b1_min,b2_deg,b2_min;

    // --------------------------------------------------------------------- //
    //                     TEST INVERSE VINCENTY ALGORITHM
    // --------------------------------------------------------------------- //
    printf("\n> TEST INVERSE VINCENTY ALGORITHM");
    printf("\n------------------------------------");
    //  Try the (inverse) vincenty algorithm; test data/output is taken from
    //+ http://www.ga.gov.au/geodesy/datums/vincenty_inverse.jsp
    //  Initialize two points on the ellipsoid.
    p1.x = -(37.0 + 57/60.0e0 + 3.72030/3600.0e0); // latitude in degrees
    p1.y = 144.0  + 25/60.0e0 + 29.52440/3600.0e0; // longtitude in degrees
    p2.x = -(37.0 + 39/60.0e0 + 10.15610/3600.0e0);// latitude in degrees
    p2.y = 143.0  + 55/60.0e0 + 35.38390/3600.0e0; // longtitude in degrees
    
    //  Transform decimal degrees to hexicondal degrees for nice printing
    decd2hexd(p1.x, a1_deg, a1_min, a1_sec);
    decd2hexd(p1.y, a2_deg, a2_min, a2_sec);
    decd2hexd(p2.x, b1_deg, b1_min, b1_sec);
    decd2hexd(p2.y, b2_deg, b2_min, b2_sec);
    printf("\nWe start with two points on the ellipsoid, with coordinates:");
    printf("\n |φ|    |%+04i %02i %10.7f|  |φ|    |%+04i %02i %10.7f|",
        a1_deg,a1_min,a1_sec,b1_deg,b1_min,b1_sec);
    printf("\n |λ|  = |%+04i %02i %10.7f|  |λ|  = |%+04i %02i %10.7f|",
        a2_deg,a2_min,a2_sec,b2_deg,b2_min,b2_sec);
    printf("\n |h|A   |        %10.7f|  |h|B   |        %10.7f|",
        0e0, 0e0);

    //  Perform the inverse Vincenty calculation, to find forward and backward
    //+ azimouths and distance between the two points.
    //  Results should be:
    //  forward azimouth:  306° 52' 5.37
    //  backward azimouth: 127° 10' 25.07
    //  distance         : 54972.271 m
    S = inverse_vincenty(deg2rad(p1.x), deg2rad(p1.y),
                         deg2rad(p2.x), deg2rad(p2.y),
                         a_for, a_bac, 1e-12);
    //  Verify results
    assert( std::abs(S-54972.271) < 1e-3 ); // millimeter accuracy
    assert( std::abs(a_for-hexd2rad(306, 52, 5.37)) < 3*1e-5 ); // .1 sec accuracy
    assert( std::abs(a_bac-hexd2rad(127, 10,25.07)) < 3*1e-5 ); // .1 sec accuracy
    
    printf("\nInverse Vincenty Formula, gives:");
    printf("\nDistance of the ellipsoidal: %+15.5f m", S);
    rad2hexd(a_for, a_deg, a_min, a_sec);
    printf("\nForward Azimouth           : %+3d %2d %8.5f", a_deg, a_min, a_sec);
    rad2hexd(a_bac, a_deg, a_min, a_sec);
    printf("\nBackward Azimouth          : %+3d %2d %8.5f", a_deg, a_min, a_sec);
    printf("\n> Everything looks OK (checked with results from");
    printf("\nhttp://www.ga.gov.au/geodesy/datums/vincenty_inverse.jsp\n");
    
    // --------------------------------------------------------------------- //
    //                     TEST VINCENTY ALGORITHM
    // --------------------------------------------------------------------- //
    //  Do the inverse from above (so results should be the same, going back
    //+ to p2)
    printf("\n> TEST DIRECT VINCENTY ALGORITHM");
    printf("\n------------------------------------");
    printf("\nWe will now perform the inverse computation, i.e. use the reults");
    printf("\nfrom the above computation and perform the direct Vincenty Formula");
    printf("\nto check results.");
    //  Compute the direct Vincenty, using point 1 and the results we got above
    new_az = direct_vincenty(deg2rad(p1.x), deg2rad(p1.y),
                             a_for, S,
                             new_lat, new_lon, 1e-12);
    // Verify results
    assert(std::abs(new_lat-deg2rad(p2.x)) < 3*1e-7 ); // .001 sec accuracy
    assert(std::abs(new_lon-deg2rad(p2.y)) < 3*1e-7 ); // .001 sec accuracy
    assert(std::abs(new_az-a_bac) < 3*1e-7 );          // .001 sec accuracy

    printf("\nFrom the direct computation, we get the results:");
    rad2hexd(new_lon, a_deg, a_min, a_sec);
    printf("\nLongtitude:  %+3d %2d %8.5f", a_deg, a_min, a_sec);
    rad2hexd(new_lat, a_deg, a_min, a_sec);
    printf("\nLatitude  :  %+3d %2d %8.5f", a_deg, a_min, a_sec);
    rad2hexd(new_az, a_deg, a_min, a_sec);
    printf("\nAzimouth  :  %+3d %2d %8.5f", a_deg, a_min, a_sec);

    printf("\nDifferences between the computed and inserted values, are:)");
    printf("\nLatitude   : %+15.5f sec.", rad2deg(std::abs(new_lat-deg2rad(p2.x)))*3600.0);
    printf("\nLongtitude : %+15.5f sec.", rad2deg(std::abs(new_lon-deg2rad(p2.y)))*3600.0);
    printf("\nAzimouth   : %+15.5f sec.", rad2deg(std::abs(new_az-a_bac))*3600.0);

    printf("\n> Everything looks OK!\n");
    
    // Get a hint of the internal accuracy of the inverse and direct transform
    double d;
    for (d = 3*1e-7; d >= 1e-15; d *= 1e-1) {
        if (  std::abs(new_lat-deg2rad(p2.x)) > d
            ||std::abs(new_lon-deg2rad(p2.y)) > d
            || std::abs(new_az-a_bac) > d ) {
            break;
        }
    }
    printf("\nInternal accuracy of the Inverse vs Direct Vincenty implementation,");
    printf("\nis %e radians, or %e arcsec\n", d, rad2deg(d)*3600e0);
    return 0;
}
