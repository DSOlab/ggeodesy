# Introduction

This is a C++ library meant to provide implementations of the most commonly used
geodetic calculations. The whole library is wrapped around the `ngpt` namespace.

Here is a list of the provided utilities:

## Reference Ellipsoids

Currently there are implementations for the (reference) ellipsoids GRS80,
WGS84 and PZ90; users can easily add more if they need to. Ellipsoid is a
class (i.e. ```ngpt::ellipsoid```) and provides access to the most common
geometric characteristics like 
```
    using namespace ngpt;
    double e2 = eccentricity_squared<ellipsoid::wgs84>();
    double b  = semi_minor<ellipsoid::grs80>();
    // normal radius of curvature at a given latitude
    double n  = N<ellipsoid::pz90>(/*latitude in radians*/);
    // meridional radii of curvature at a given latitude
    double m  = M<ellipsoid::pz90>(/*latitude in radians*/);
```
