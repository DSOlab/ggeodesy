#include "meridian_arc_length.hpp"
#ifdef DEBUG
#include <iostream>
#endif

/// @brief Meridian arc length on ellipsoid from the equator to latitude
/// a: semi-major
/// Implementation using a binomial series with respect to e^2 obtained by
/// truncating the expansion at order e^10
/// K. Kawase, "A General Formula for Calculating Meridian Arc Length and its
/// Application to Coordinate Conversion in the Gauss-Krüger Projection"
double ngpt::core::meridian_arc_length_impl1(double a, double f,
                                             double lat) noexcept {
  const double e2 = eccentricity_squared(f);
  const double e4(e2 * e2), e6(e4 * e2), e8(e6 * e2), e10(e8 * e2);
  const double fac = a * (1e0 - e2);
  const double C1 =
      1e0 +
      e2 * (3e0 / 4 + e2 * (45e0 / 64 + e2 * (175e0 / 256 +
                                              e2 * (11025e0 / 16384 +
                                                    (43659e0 / 65536) * e2))));
  const double C2 =
      e2 *
      (3e0 / 4 +
       e2 * (15e0 / 16 + e2 * (525e0 / 512 +
                               e2 * (2205e0 / 2048 + (72765e0 / 65536) * e2))));
  const double C3 =
      e4 * (15e0 / 64 +
            e2 * (105e0 / 256 + e2 * (2205e0 / 4096 + (10395e0 / 16384) * e2)));
  const double C4 =
      e6 * (35e0 / 512 + e2 * (315e0 / 2048 + (31185e0 / 131072) * e2));
  const double C5 = e8 * (315e0 / 16384 + (3465e0 / 65536) * e2);
  const double C6 = e10 * (693e0 / 131072);
  const double sin2f(std::sin(2e0 * lat)), sin4f(std::sin(4e0 * lat)),
      sin6f(std::sin(6e0 * lat)), sin8f(std::sin(8e0 * lat)),
      sin10f(std::sin(10e0 * lat));
  return fac * (C1 * lat - C2 * sin2f / 2e0 + C3 * sin4f / 4e0 -
                C4 * sin6f / 6e0 + C5 * sin8f / 8e0 - C6 * sin10f / 10e0);
}

/// @brief Meridian arc length on ellipsoid from the equator to latitude
/// a: semi-major
/// This is the Bessel’s formula implementation
/// K. Kawase, "A General Formula for Calculating Meridian Arc Length and its
/// Application to Coordinate Conversion in the Gauss-Krüger Projection"
double ngpt::core::meridian_arc_length_impl2(double a, double f,
                                             double lat) noexcept {
  const double b = ngpt::core::semi_minor(a, f);
  // third flattening
  const double n = (a - b) / (a + b);
  const double n2(n * n), n3(n2 * n);
  const double onemn = (1e0 - n);
  const double onepn = (1e0 + n);
  const double a1mn1pn = a * onemn * onepn;
  const double C1 = 1e0 + n2 * (9e0 / 4 + (225e0 / 64) * n2);
  const double C2 = n * (1e0 + n2 * (15e0 / 8 + (175e0 / 64) * n2));
  const double C3 = n2 * (1e0 + (7e0 / 4) * n2);
  const double C4 = n3 * (1e0 + (27e0 / 16) * n2);
  const double sin2f(std::sin(2e0 * lat)), sin4f(std::sin(4e0 * lat)),
      sin6f(std::sin(6e0 * lat));
  double L = C1 * lat - (3e0 * sin2f * C2) / 2e0 + (15e0 * C3 * sin4f) / 16e0 -
             (35e0 * C4 * sin6f) / 48e0;
  L *= a1mn1pn;
  L *= onemn;
  return L;
}

/// @brief Meridian arc length on ellipsoid from the equator to latitude
/// a: semi-major
/// This is the Helmert's formula implementation
/// K. Kawase, "A General Formula for Calculating Meridian Arc Length and its
/// Application to Coordinate Conversion in the Gauss-Krüger Projection"
double ngpt::core::meridian_arc_length_impl3(double a, double f,
                                             double lat) noexcept {
  const double b = ngpt::core::semi_minor(a, f);
  // third flattening
  const double n = (a - b) / (a + b);
  const double n2(n * n), n3(n2 * n), n4(n3 * n);
  const double C1 = 1e0 + n2 * (1e0 / 4e0 + n2 * (1e0 / 64e0));
  const double C2 = n * (3e0 / 2e0 - n2 * (3e0 / 16e0));
  const double C3 = n2 * (15e0 / 16e0 - n2 * (15e0 / 64e0));
  const double C4 = 35e0 * n3 / 48e0;
  const double C5 = 315e0 * n4 / 512e0;
  const double sin2f(std::sin(2e0 * lat)), sin4f(std::sin(4e0 * lat)),
      sin6f(std::sin(6e0 * lat)), sin8f(std::sin(8e0 * lat));
  double L = a * (C1 * lat - C2 * sin2f + C3 * sin4f - C4 * sin6f + C5 * sin8f);
  return L / (1e0 + n);
}
