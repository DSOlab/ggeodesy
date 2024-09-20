# C++ Geodetic Library

[![clang-format Check](https://github.com/DSOlab/ggeodesy/actions/workflows/clang-format-check.yml/badge.svg)](https://github.com/DSOlab/ggeodesy/actions/workflows/clang-format-check.yml)
[![Linux CI build](https://github.com/DSOlab/ggeodesy/actions/workflows/cpp-linux-build.yml/badge.svg)](https://github.com/DSOlab/ggeodesy/actions/workflows/cpp-linux-build.yml)

# Introduction

This is a C++ library meant to provide implementations of the most commonly used
geodetic calculations. The whole library is wrapped around the `dso` namespace.

# Dependancies 

This library uses the [eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) 
library for basic matrix manipulation and linear algebra.

# Compilation / Installation

Building the library requires [cmake](https://cmake.org/)
Supposing you are located in the top-level directory:

```
## to build in a folder named "build":
$> cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
$> cmake --build build --target all --config Release -- -j4
## (Optional) run tests
$> ctest --test-dir build
## Install, system-wide (needs root)
$> cd build && sudo make install
```

# The Library

## Coordinate Types

A point $P$ can be defined can be described (in 3-D space) by any of the following 
coordinate types:

 * **Cartesian**, using $P = (x,y,z)$; unless otherwise stated, in SI units (i.e. meters)
 * **Geodetic**, using $P= (\lambda , \phi , h)$, $\lambda$ denoting the 
  longitude in range $-\pi \le \lambda \le \pi$, $\phi$ denoting the geodetic
  latitude in range $\frac{\pi}{2} \le \phi \le \frac{\pi}{2}$ and $h$ denotes the 
  ellipsoidal height. Unless otherwise stated, ellipsoidal coordinate sets are 
  given/derived in **this order** (i.e. $(\lambda, \phi, h)$) in units of radians 
  and meters. Note that **geodetic** coordinates are based on a reference ellipsoid.
 * **Spherical**

## Coordinate Transformations

- Test : Geodetic -> Cartesian -> Geodetic results in max discrepancies (between
  the input and ouput geodetic coordinates) in the range:
   $max\delta \phi \approx 1e^{-10} arcsec$, $max\delta \lambda \approx 5e^{-11} arcsec$ 
   and $max\delta height \approx 4e^{-9} m$. See [here](test/unit/geodetic.cpp)).

- Test : Spherical -> Cartesian -> Spherical results in max discrepancies (between
  the input and ouput spherical coordinates) in the range:
   $max\delta \phi _{geocentric} \approx 1e^{-8} arcsec$, $max\delta \lambda \approx 5e^{-11} arcsec$ 
   and $max\delta height \approx 2e^{-9} m$. See [here](test/unit/spherical.cpp)).
