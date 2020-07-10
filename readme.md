# C++ Geodetic Library

[![Build Status](https://travis-ci.com/xanthospap/ggeodesy.svg?branch=master)](https://travis-ci.com/xanthospap/ggeodesy)

## Introduction

This is a C++ library meant to provide implementations of the most commonly used
geodetic calculations. The whole library is wrapped around the `ngpt` namespace.

## Compilation / Installation

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
```


## Verify & Test

In the `ggeodesy/test` folder you should be able to see a list of executables; run
`ggeodesy/test/testGeodesy` to validate the library.

After a succesefull installation, users should have:

1. all library header files in `/usr/local/include/ggeodesy/`
2. the library (both a static and shared) in `/usr/local/lib/`

Link, include and have fun!

## The Library

Here is a list of the provided utilities:

### Reference Ellipsoids

Currently there are implementations for the (reference) ellipsoids:

  - [GRS80](https://en.wikipedia.org/wiki/Geodetic_Reference_System_1980),
  - [WGS84](https://en.wikipedia.org/wiki/World_Geodetic_System) and 
  - [PZ90](https://eng.mil.ru/files/PZ-90.11_final-v8.pdf)

Reference ellipsoids can be used in either on of two ways:

 - via using the `enum` class `ngpt::ellipsoid` (e.g. `ngpt::ellipsoid::grs80`), or
 - via the class `ngpt::Ellipsoid`

Users can easily add more reference ellipsoids if they need to, via constructing 
an `Ellipsoid` instance with the wanted parameters (aka semi-major axis and flattening), e.g.

```
    using namespace ngpt;
    Ellipsoid myell = Ellipsoid(6378136.0e0/*semi-major axis*/, 1/298.25784/*flattening factor*/);
    // use the created ellipsoid in some way ....
    double semi_major = myell.semi_major();
    double N = myell.N(some_latitude);
    // ....
```

Users can also expand the source code to add a reference ellipsoid in the `enum` (class)
`ellipsoid`. In this case, you will need one entry in the `ellipsoid` enum and a
respective specialization of the `template <ellipsoid E> struct ellipsoid_traits {};` 
class.

Note that most of the constructors and function (for the `Ellipsoid` class and the 
`ellipsoid` enum are `constexpr`. Hence, the following code is computed at compile-time:

```
  using namespace ngpt;
  constexpr auto wgs84 = Ellipsoid(ellipsoid::wgs84);
  constexpr auto grs80 = Ellipsoid(ellipsoid_traits<ellipsoid::grs80>::a,
                                   ellipsoid_traits<ellipsoid::grs80>::f);
  constexpr auto pz90  = Ellipsoid(ellipsoid::pz90);

  static_assert(wgs84.eccentricity_squared() == 
                eccentricity_squared<ellipsoid::wgs84>());
  static_assert(grs80.semi_minor() == semi_minor<ellipsoid::grs80>());
  static_assert(std::abs(pz90.eccentricity_squared()-0.0066943662)<1e-9);
```

For more information on how to use the reference ellipsoids, see e.g. [test_ellipsoid.hpp](https://github.com/xanthospap/ggeodesy/blob/master/test/test_ellipsoid.cpp).

### Coordinate Transformations

The following coordinate transformations are provided (for points on some reference ellipsoid):
* Cartesian to Ellipsoidal (aka [x, y, z] to [φ, λ, height])
* Ellipsoidal to Cartesian (aka [φ, λ, height] to [x, y, z])
* Cartesian to Topocentric (aka [δx, δy, δz] to [north, east, up])
* Topocentric to Cartesian (aka  [north, east, up] to [δx, δy, δz])

### How to use the library (TODO)

### Namespaces

The whole of the library is wrapped around the `ngpt` namespace

### Linking

- static
- dynamic

## Documentation & Library API (TODO)

- build dox with doxygen (or link to dox)

## TODO

## Bugs & Maintanance
Xanthos, xanthos@mail.ntua.gr
Mitsos, danast@mail.ntua.gr
