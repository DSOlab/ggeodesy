# C++ Geodetic Library

[![Build Status](https://travis-ci.com/xanthospap/ggeodesy.svg?branch=master)](https://travis-ci.com/xanthospap/ggeodesy)

## Introduction

This is a C++ library meant to provide implementations of the most commonly used
geodetic calculations. The whole library is wrapped around the `ngpt` namespace.

## Compilation / Installation

Source code is ISO C++17. Compilation should be trivial using any C++ compiler
[supporting the c++17](https://en.wikipedia.org/wiki/C%2B%2B17#Compiler_support) 
standard (option `-std=c++17` in gcc and clang).

> This software is meant to be implemented on Unix-type OS's. No effort will be
> undertaken for compatibility with other OS types.

To compile the library, just follow the basic steps: (*note that the library is still at development phase so users need to configure the project before compiling*)

For the following, `ROOTDIR` will be the root directory of this repository,
aka the directory under which `/src`, `/test` and `/doc` folders live.

~~**If you do not need the DEBUG version** (which most probably you don't), create the `Makefile.am` templates. This means that you
should rename [Makefile.am.production](src/Makefile.am.production) and [Makefile.am.production](test/Makefile.am.production) to
`src/Makefile.am` and `test/Makefile.am` respectively, that is:~~

To prepare the required files for compilation (that is the `Makefile.am` in each 
of the relevant folders) you need to run the script [install_setup.py](install_setup.py). 
You can use the `-h` switch to see the help message, but in most cases the 
command `./install_setup.py -c production` will suffice.

If needed (that is you are not running the script from `ROOTDIR`) specify the 
`ROOTDIR` path via the `-d` switch.

Then run Autotools and compile:

```bash
autoreconf -if
./configure
make
sudo make install
```

## Verify & Test

~~In the `ggeodesy/test` folder you should be able to see a list of executables; run `ggeodesy/test/testGeodesy` to validate the library.~~

After a succesefull installation, users should have:

1. all library header files in `/usr/local/include/ggeodesy/`
2. the library (both a static and shared) in `/usr/local/lib/`

To run a validity check, just run: `make check` at the root directory. Hopefully, 
you 'll see all checks passing!

Link, include and have fun!

## The Library

Here is a list of the provided utilities:

### Radians, Degrees and Relevant Trnasformations

The library includes functions for transforming between radians, decimal degrees 
and hexicondal degrees. Namely:

  - `ngpt::deg2rad` transforms radians to decimal degrees
  - `ngpt::decd2hexd` transforms decimal degrees to hexicondal degrees
  - `ngpt::rad2deg` transforms decimal degrees to radians
  - `ngpt::rad2hexd` transforms radians to hexicondal degrees
  - `ngpt::hexd2decd` transforms hexicondal degrees to decimal degrees
  - `ngpt::hexd2rad` transforms hexicondal degrees to radians

Additionaly the function `ngpt::normalize_angle` can normalize a given angle 
to a specified range (e.g. in range -π to π).

All of the above functions are defined in the header file [geodesy.hpp](https://github.com/xanthospap/ggeodesy/blob/master/src/geodesy.hpp). For usage examples, see 
[test_units.hpp](https://github.com/xanthospap/ggeodesy/blob/master/test/test_units.cpp) and [test_angle_normalization.hpp](https://github.com/xanthospap/ggeodesy/blob/master/test/test_angle_normalization.cc).

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

```cpp
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

```cpp
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

### Meridian Arc Length

As is well known in geodesy, the meridian arc length S(\f$\phi\f$) on the earth ellipsoid 
from the equator to the geographic latitude \f$\phi\f$ includes an elliptic integral and 
cannot be expressed explicitly using a combination of elementary functions. The library 
provides three implementations for computing the meridian arc length using approximations.

* `ngpt::core::meridian_arc_length_impl1` : implementation using a binomial series 
with respect to e^2 obtained by truncating the expansion at order e^10

* `ngpt::core::meridian_arc_length_impl2` : the Bessel’s formula implementation

* `ngpt::core::meridian_arc_length_impl3` : the Helmert's formula implementation

The above implementations agree within the following precision limits:

* impl1 vs impl2 : 3e-5 meters
* impl2 vs impl3 : 1e-6 meters
* impl2 vs impl3 : 3e-5 meters

Users can select one of the algorithms via the (optional) intput parameter `alg` 
(using values in range [0,2]) when calling the function 
`template<ellipsoid E> double meridian_arc_length(double lat, int alg=0)` (aka 
by default the function will use the 

More information and implementation details can be found in [Kawase, 2011](#kawase). 
Usage examples of the algorithm(s) can be found in the file 
[test_meridian_arc.hpp](https://github.com/xanthospap/ggeodesy/blob/master/test/test_meridian_arc.cpp)

Note that if we only want the __arc length of an infinitesimal element of the meridian__ the 
computation is way more straight-forward [Meridian Arc](#meridian_arc_wiki); for this computation users may use the 
function (template) `infinitesimal_meridian_arc`.

### How to use the library (TODO)

### Namespaces

The whole of the library is wrapped around the `ngpt` namespace

### Linking

- static
- dynamic

## Documentation & Library API (TODO)

- build dox with doxygen (or link to dox)

## FAQ

## TODO

## Bugs & Maintanance
Xanthos, xanthos@mail.ntua.gr
Mitsos, danast@mail.ntua.gr

## References 

<a id="kawase"></a> K. Kawase, A General Formula for Calculating Meridian Arc Length and its Application to Coordinate 
Conversion in the Gauss-Krüger Projection, Bulletin of the Geospatial Information Authority of Japan, Vol.59 December, 2011

<a id="meridian_arc_wiki"></a> Wiki page on [Meridian Arc](https://en.wikipedia.org/wiki/Meridian_arc)

<a id="moritz-grs80"></a>H. Moritz, GEODETIC REFERENCE SYSTEM 1980, available
[here](https://geodesy.geology.ohio-state.edu/course/refpapers/00740128.pdf)
