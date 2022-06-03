# C++ Geodetic Library

[![Build Status](https://travis-ci.com/xanthospap/ggeodesy.svg?branch=master)](https://travis-ci.com/xanthospap/ggeodesy)

## Introduction

This is a C++ library meant to provide implementations of the most commonly used
geodetic calculations. The whole library is wrapped around the `dso` namespace.

## Compilation / Installation

> Since December 2021, the build system has been changed from 
> [GNU Autotools](https://www.gnu.org/software/automake/manual/html_node/Autotools-Introduction.html) 
> to [scons](https://scons.org/). 


### TL;DR

First off, clone the project.

To buid the project just use [scons](https://scons.org/):
```
scons [OPTIONS]
```
in the top-level directory. The optional `[OPTIONS]` argument, can be any/multiple of:
 * `debug=1` to trigger a `DEBUG` build, 
 * `boost=1` to trigger a build including [boost geometry](https://www.boost.org/doc/libs/1_74_0/libs/geometry/doc/html/index.html) 
    comparisson/test programs
 * `eigen=1` to use [eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) 
    library for vector/matrix operations. **Note that because the library uses extensively 
    template code, if you include `ggeodesy` header files in another project you should 
    also include the `-DUSE_EIGEN` compilation flag when building the dependent project**


Installation is trivial, use:
```
sudo scons install
```

### testing against boost/geometry

The folder [boost](boost) includes source code for testing the algorithms in 
ggeodesy against the ones implemented in 
(boost/geometry)[https://www.boost.org/doc/libs/1_74_0/libs/geometry/doc/html/index.html].
To compile these, you will need the respective developement files. Installing 
boost is usually pretty trivial; relevant documentation can be found on the 
(boost webpage)[https://www.boost.org/].

Once you have downloaded `boost/geometry` you can include the [boost](boost) 
source code in the build process via `scons boost=1`


## Verify & Test

TODO

## Accuracy/Precision and Floating Point Numbers Considerations

For testing purposes, we need to compare (floating point) results, aka check if 
two floating point numbers are __equal__. But what does equal mean? Well, in 
this case, for most of the test programs we define the floating point 
comparisson (when comparing numbers other than zero) as:

```
diff = abs(a - b)
diff < max(abs(a), abs(b)) * tolerance or diff < tolerance
```

where __tolerance__ is defined as: `std::numeric_limits<double>::epsilon()`. 
That is, we use a tolerance value proportional to max(a,b) when comparing two 
(floating point) numbers.
(see also [this SO thread](https://stackoverflow.com/questions/17333/what-is-the-most-effective-way-for-float-and-double-comparison), 
and [this also](https://stackoverflow.com/questions/48133572/what-can-stdnumeric-limitsdoubleepsilon-be-used-for)).

This comparisson function is defined in the 
[test_help.hpp](https://github.com/xanthospap/ggeodesy/blob/master/test/test_help.hpp) 
file and used throughout the test/check programs.

## The Library

Here is a list of the provided utilities:

### General Notes

In general, when refering to coordinate arrays/vector, the elements are in the 
following order:

  * for cartesian vector `[x, y, z]`
  * for ellipsoidal/geodetic (longitude, latitude, height) `[longitude, latitude, height]`
  * topocentric `[east, north, up]`

### Radians, Degrees and Relevant Transformations

The library includes functions for transforming between radians, decimal degrees 
and hexicondal degrees. Namely:

  - `ngpt::deg2rad` transforms radians to decimal degrees
  - `ngpt::decd2hexd` transforms decimal degrees to hexicondal degrees
  - `ngpt::rad2deg` transforms decimal degrees to radians
  - `ngpt::rad2hexd` transforms radians to hexicondal degrees
  - `ngpt::hexd2decd` transforms hexicondal degrees to decimal degrees
  - `ngpt::hexd2rad` transforms hexicondal degrees to radians
  - `ngpt::rad2sec` convert radians to seconds (of degree)

Additionaly the function `ngpt::normalize_angle` can normalize a given angle 
to a specified range (e.g. in range -π to π).

All of the above functions are defined in the header file [units.hpp](https://github.com/xanthospap/ggeodesy/blob/master/src/units.hpp). For usage examples, see 
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

The whole of the library is wrapped around the `dso` namespace

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
