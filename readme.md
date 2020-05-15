# Introduction

This is a C++ library meant to provide implementations of the most commonly used
geodetic calculations. The whole library is wrapped around the `ngpt` namespace.

# Compilation / Installation

Source code is ISO C++17. Compilation should be trivial using any gcc version 
supporting the c++17 standard (option `-std=c++17`).

> This software is meant to be implemented on Unix-type OS's. No effort will be
> undertaken for compatibility with other OS types.

To compile the library, just follow the basic steps: (*note that the library is still at development phase so users need to configure the project before compiling*)

**If you do not need the DEBUG version** (which most probably you don't), create the `Makefile.am` templates. This means that you
should rename [Makefile.am.production](src/Makefile.am.production) and [Makefile.am.production](test/Makefile.am.production) to
`src/Makefile.am` and `test/Makefile.am` respectively, that is:
```
mv src/Makefile.am.production src/Makefile.am
mv test/Makefile.am.production test/Makefile.am
```

Then run Autotools and compile:

```
autoreconf -if
./configure
make
make install


## Verify & Test

In the `ggeodesy/test` folder you should be able to see a list of executables; run
`ggeodesy/test/testGeodesy` to validate the library.

After a succesefull installation, users should have:

1. all library header files in `/usr/local/include/ggeodesy/`
2. the library (both a static and shared) in `/usr/local/lib/`

Link, include and have fun!

# The Library

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

## How to use the library (TODO)

### Namespaces

- namespace `iers2010`
- namespace `iers2010::dhtide`
- namespace `iers2010::hisp`
- namespace `iers2010::oeop`

### Linking

- static
- dynamic

## Documentation & Library API (TODO)

- build dox with doxygen (or link to dox)

## TODO

## Bugs & Maintanance
Xanthos, xanthos@mail.ntua.gr
Mitsos, danast@mail.ntua.gr
