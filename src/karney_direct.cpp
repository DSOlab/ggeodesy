#include "vincenty.hpp"
#include "geoconst.hpp"
#include "geodesy.hpp"
#include <cmath>
#include <stdexcept>
#ifdef DEBUG
#include "units.hpp"
#include <fenv.h>
#endif

using ngpt::D2PI;
using ngpt::DPI;

/// @brief Direct Karney formula for the direct geodesic problem.
///
/// Given an initial point (lat1, lon1), an initial azimuth, a1, and a distance,
/// s, along the geodesic the problem is to find the end point (lat2, lon2)
/// and azimuth, a2.
///
/// @param[in]  lat1        Latitude of point 1 (radians)
/// @param[in]  lon1        Longtitude of point 1 (radians)
/// @param[in]  a1          Azimouth from point 1 to point 2 (radians)
/// @param[in]  s           Ellipsoid distance (along the geodesic) from point
///                         1 to point 2 (meters)
/// @param[in]  semi_major  Semi-major axis of reference ellipsoid (meters)
/// @param[in]  semi_minor  Semi-minor axis of reference ellipsoid (meters)
/// @param[in]  flattening  Flattening of reference ellipsoid
/// @param[out] lat2        Latitude of point 2 (radians)
/// @param[out] lon2        Longtitude of point 2 (radians)
/// @param[in]  convergence_limit Convergence limit (radians). 10e-12
///                         corresponds to approximately 0.06mm
/// @return     Azimouth from point 2 to point 1 (radians)
///
/// @see https://link.springer.com/article/10.1007/s00190-012-0578-z
///      Karney 2012, Algorithms for geodesics, Journal of Geodesy volume 87, 
///      pages43–55(2013)
double ngpt::core::direct_karney(double lat1, double lon1, double a1,
                                   double s, double semi_major,
                                   double flattening, double semi_minor,
                                   double &lat2, double &lon2) {

  const double a = semi_major;
  const double f = flattening;
  const double b = semi_minor;

  printf("\nKarney");
  // Solve triangle NEA
  double beta = ngpt::core::reduced_latitude(lat1, f);
  printf("\nbeta=%15.10f", rad2deg(beta));
  const double sinbeta = std::sin(beta);
  const double cosbeta = std::cos(beta);
  const double sinalpha = std::sin(a1);
  const double cosalpha = std::cos(a1);
  double cmplx_norm = std::sqrt(cosalpha*cosalpha + (sinalpha*sinbeta)*(sinalpha*sinbeta));
  const double alpha0 = std::atan2(sinalpha*cosbeta, cmplx_norm);
  printf("\nalpha0=%15.10f", rad2deg(alpha0));
  const double sinalpha0 = std::sin(alpha0);
  //const double cosalpha0 = std::cos(alpha0);
  const double sigma = std::atan2(sinbeta, cosalpha*cosbeta);
  printf("\nsigma=%15.10f", rad2deg(sigma));
  const double sinsigma = std::sin(sigma);
  const double cossigma = std::cos(sigma);
  const double omega = std::atan2(sinalpha0*sinsigma, cossigma);
  printf("\nomega=%15.10f", rad2deg(omega));
  // cmplx_norm = std::sqrt((cosalpha0*cossigma)*(cosalpha0*cossigma)+sinalpha0*sinalpha0);
  // beta = std::atan2(cosalpha0*sinsigma, cmplx_norm);
  // alpha = std::atan2(sinalpha0, cosalpha0*cossigma);

  lat2=lon1;
  lon2=a;
  lon1=b;
  return s;
}
