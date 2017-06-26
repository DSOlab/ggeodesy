#include <stdio.h>
#include <cassert>
#include "geodesy.hpp"
#include "car2top.hpp"
#include "car2ell.hpp"

struct Point { double x,y,z; };

// distance between two points
double
distance(const Point& p1, const Point& p2) noexcept
{ return std::sqrt( (p1.x-p2.x)*(p1.x-p2.x)
                  + (p1.y-p2.y)*(p1.y-p2.y)
                  + (p1.z-p2.z)*(p1.z-p2.z) );
}

// norm of vector
double
distance(const Point& p1) noexcept
{ return std::sqrt( p1.x*p1.x
                  + p1.y*p1.y
                  + p1.z*p1.z );
}

// precision for comparing two doubles for equality
constexpr double prec {1e-10};

using namespace ngpt;

int main ()
{
    /// a (reference) point on ellipsoid.
    Point p1 {4595220.02261, 2039434.13622, 3912625.96044};
    Point p2 {4595220.02261+5.12345,
              2039434.13622+5.12345,
              3912625.96044+5.12345};
    
    // --------------------------------------------------------------------- //
    //                     TEST COORDINATE TRANSFORMATIONS
    // --------------------------------------------------------------------- //
    printf("\n> TEST COORDINATE TRANSFORMATIONS\n");
    printf("------------------------------------\n");
    printf("\n\nWe start with the following two points (cartesian):");
    printf("\n| x |   |%+15.5f| | x |   |%+15.5f|", p1.x, p2.x);
    printf("\n| y | = |%+15.5f| | y | = |%+15.5f|", p1.y, p2.y);
    printf("\n| z |A  |%+15.5f| | z |B  |%+15.5f|", p1.z, p2.z);
    
    printf("\n\nCartesian and Topocentric vectors, between pA and pB are:");
    Point p3 {p1.x-p2.x, p1.y-p2.y, p1.z-p2.z}; // [dx,dy,dz] vector
    Point p4; // will hold the topocentric vector, [north, east, up]
    ngpt::car2top<ellipsoid::grs80>(p1.x, p1.y, p1.z, p2.x, p2.y, p2.z, p4.x, p4.y, p4.z);
    printf("\n| δx |   |%+15.10f|  | n |   |%+15.10f|", p3.x, p4.x);
    printf("\n| δy | = |%+15.10f|  | e | = |%+15.10f|", p3.y, p4.y);
    printf("\n| δz |   |%+15.10f|  | u |   |%+15.10f|", p3.z, p4.z);

    double d1,d2,d3;
    d1 = distance(p1, p2);
    d2 = distance(p3);
    d3 = distance(p4);
    printf("\n\nDistance between points p1 and p2, calculated from cartesian");
    printf("\ncomponents, is                        : %15.7f", d1);
    printf("\nNorm of the cartesian diff. vector is : %15.7f", d2);
    printf("\nNorm of the topocentric vector is     : %15.7f", d3);
    assert( std::abs(d1-d2) < prec && std::abs(d2-d3) < prec );

    printf("\n\nFrom the topocentric [N, E, U] vector, we can calculate the");
    printf("\ndistance, azimouth and zenith angle between the two points:");
    double dist, azi, zen;
    top2daz(p4.x, p4.y, p4.z, dist, azi, zen);
    printf("\ndistance   = %15.7f (m)", dist);
    printf("\nazimouth   = %15.10f deg.", rad2deg(azi) );
    printf("\nzenith ang.= %15.10f deg.", rad2deg(zen) );
    assert( std::abs(dist-d1) < prec );

    printf("\n\nFrom the topocentric vector we can go back to cartesian, if we");
    printf("\nknow the coordinates (actually latitude and longtitude of the");
    printf("\nreference point. Hence, actually compute the cartesian components");
    printf("\nof the second, ending point");
    Point p5, // ellipsoidal coordinates of p1
          p6; // cartesian coordinates of p2
    car2ell<ellipsoid::grs80>(p1.x, p1.y, p1.z, p5.x, p5.y, p5.z);
    top2car(p4.x, p4.y, p4.z, p5.x, p5.y, p6.x, p6.y, p6.z);
    printf("\n| φ |   |%+15.10f| | δx |    |%+15.5f|", rad2deg(p5.x), p6.x);
    printf("\n| λ | = |%+15.10f| | δy | =  |%+15.5f|", rad2deg(p5.y), p6.y);
    printf("\n| h |A  |%+14.5f | | δz |AB  |%+15.5f|", p5.z, p6.z);
    printf("%15.10f",std::abs(p6.x+p3.x));
    printf("%15.10f",std::abs(p6.y+p3.y));
    printf("%15.10f",std::abs(p6.z+p3.z));
    assert(  std::abs(p6.x+p3.x) < prec
          && std::abs(p6.y+p3.y) < prec
          && std::abs(p6.z+p3.z) < prec );

    printf("\n");
    return 0;
}
