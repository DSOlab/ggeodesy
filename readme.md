# Introduction

This is a C++ library meant to provide implementations of the most commonly used
geodetic calculations. The whole library is wrapped around the `ngpt` namespace.

# Installation

* Clone the repository
`git clone https://xanthos@bitbucket.org/xanthos/ggeodesy.git` into a local dir 
(let's call it ggeodesy)

* Setup the Makefile.am files in the src and test directories
```
cd ggeodesy
mv src/Makefile.am.production src/Makefile.am
mv test/Makefile.am.production test/Makefile.am
```

* Configure and Compile the software
```
autoreconf -if
./configure
make
sudo make install
```

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

## Coordinate Transformations

The following coordinate transformations are provided (for points on some reference ellipsoid):
* Cartesian to Ellipsoidal (aka [x, y, z] to [φ, λ, height])
* Ellipsoidal to Cartesian (aka [φ, λ, height] to [x, y, z])
* Cartesian to Topocentric (aka [δx, δy, δz] to [north, east, up])
* Topocentric to Cartesian (aka  [north, east, up] to [δx, δy, δz])

## Other thing ... (todo)
