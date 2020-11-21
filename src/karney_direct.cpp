#include "geoconst.hpp"
#include "geodesy.hpp"
#include "vincenty.hpp"
#include <cmath>
#include <stdexcept>
#ifdef DEBUG
#include "units.hpp"
#include <fenv.h>
#endif

using ngpt::D2PI;
using ngpt::DPI;

// relative floating point comparisson;
// see https://en.cppreference.com/w/cpp/types/numeric_limits/epsilon
bool double_almost_equal(double x, double y, int ulp=2)
{
  // the machine epsilon has to be scaled to the magnitude of the values used
  // and multiplied by the desired precision in ULPs (units in the last place)
  return std::fabs(x-y) <= std::numeric_limits<double>::epsilon() * std::fabs(x+y) * ulp
    // unless the result is subnormal
    || std::fabs(x-y) < std::numeric_limits<double>::min();
}

// distance integral expanded in a Fourier series, Karney eq. 15
// note that this is only the part in parenthesis, aka needs to be
// multiplied by A1
double series_expansion_I(double epsilon, double sigma) noexcept {
  const double e = epsilon;
  const double e2 = epsilon * epsilon;
  const double e3 = e2 * epsilon;
  const double e4 = e3 * epsilon;
  const double e5 = e4 * epsilon;
  const double C[6] = {
      e * (-0.5e0 + e2 * ((3e0 / 16e0) - (1e0 / 32e0) * e2)),
      e2 * (-(1e0 / 16e0) + e2 * ((1e0 / 32e0) - (9e0 / 2048e0) * e2)),
      e3 * (-(1e0 / 48e0) + (3e0 / 256e0) * e2),
      e4 * (-(5e0 / 512e0) + (3e0 / 512e0) * e2),
      -(7e0 / 1280e0) * e5,
      -(7e0 / 2048e0) * e5 * e};
  double sum = 0e0;
  for (int l = 0; l < 6; l++)
    sum += C[l] * std::sin(2e0 * static_cast<double>(l + 1) * sigma);
  return sum + sigma;
}

// determine σ given s via Lagrange reversion, Karney eq. 20
double series_expansion_S(double epsilon, double tau) noexcept {
  const double e = epsilon;
  const double e2 = epsilon * epsilon;
  const double e3 = e2 * epsilon;
  const double e4 = e3 * epsilon;
  const double e5 = e4 * epsilon;
  const double C[6] = {
      e * (0.5e0 + e2 * (-(9e0 / 32e0) + (205e0 / 1536e0) * e2)),
      e2 * ((5e0 / 16e0) + e2 * (-(37e0 / 96e0) + (1335e0 / 4096e0) * e2)),
      e3 * ((29e0 / 96e0) - (75e0 / 128e0) * e2),
      e4 * ((539e0 / 1536e0) - (2391e0 / 2560e0) * e2),
      e5 * (3467e0 / 7680e0),
      e5 * (38081e0 / 61440e0) * e};
  double sum = 0e0;
  for (int l = 0; l < 6; l++)
    sum += C[l] * std::sin(2e0 * static_cast<double>(l + 1) * tau);
  return sum + tau;
}

// Fourier series of longtitude equation, Karney eq. 23
// Note: we need to multiply this with A3 to get I3(σ)
// Note: we compute the C3j factors first; then compute the sum in eq. 23 for
//       two independent sigma values given in the *sigmas vector
void series_expansion_I3(double epsilon, double n, const double* sigmas, double* sums) noexcept {
  const double e = epsilon;
  const double e2 = epsilon * epsilon;
  const double e3 = e2 * epsilon;
  const double e4 = e3 * epsilon;
  const double e5 = e4 * epsilon;
  const double n2 = n * n;
  const double C[6] = {
      e * ((0.25e0 - 0.25e0 * n) +
           e * (((1e0 / 8e0) - (1e0 / 8e0) * n2) +
                e * (((3e0 / 64e0) + (3e0 / 64e0) * n - (1e0 / 64e0) * n2) +
                     e * (((5e0 / 128e0) + (1e0 / 64e0) * n) +
                          e * (3e0 / 128e0))))),
      e2 *
          (((1e0 / 16e0) - (3e0 / 32e0) * n + (1e0 / 32e0) * n2) +
           e * (((3e0 / 64e0) - (1e0 / 32e0) * n - (3e0 / 64e0) * n2) +
                e * (((3e0 / 128e0) + (1e0 / 128e0) * n) + e * (5e0 / 256e0)))),
      e3 * (((5e0 / 192e0) - (3e0 / 64e0) * n + (5e0 / 192e0) * n2) +
            e * (((3e0 / 128e0) - (5e0 / 192e0) * n) + e * (7e0 / 512e0))),
      e4 * (((7e0 / 512e0) - (7e0 / 256e0) * n) + e * (7e0 / 512e0)),
      e5 * (21e0 / 2560e0)};
  double sum = 0e0;
  for (int l = 0; l < 6; l++)
    sum += C[l] * std::sin(2e0 * static_cast<double>(l + 1) * sigmas[0]);
  sums[0] = sum + sigmas[0];
  sum =0e0;
  for (int l = 0; l < 6; l++)
    sum += C[l] * std::sin(2e0 * static_cast<double>(l + 1) * sigmas[1]);
  sums[1] = sum + sigmas[1];
  return;
}

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
double ngpt::core::direct_karney(double lat1, double lon1, double a1, double s,
                                 double semi_major, double flattening,
                                 double semi_minor, double &lat2,
                                 double &lon2) {

  const double a = semi_major;
  const double f = flattening;
  const double b = semi_minor;

  //printf("\nKarney");
  // Solve triangle NEA
  double beta1 = ngpt::core::reduced_latitude(lat1, f);
  //printf("\nbeta1=%15.11f", rad2deg(beta1));
  const double sinbeta1 = std::sin(beta1);
  const double cosbeta1 = std::cos(beta1);
  const double sinalpha = std::sin(a1);
  const double cosalpha = std::cos(a1);
  double cmplx_norm = std::sqrt(cosalpha * cosalpha +
                                (sinalpha * sinbeta1) * (sinalpha * sinbeta1));
  const double alpha0 = std::atan2(sinalpha * cosbeta1, cmplx_norm);
  //printf("\nalpha0=%15.11f", rad2deg(alpha0));
  const double sinalpha0 = std::sin(alpha0);
  double __sigma1 = std::atan2(sinbeta1, cosalpha * cosbeta1);
  // For equatorial geodesics (φ1= 0 and α1=π/2), Eq. (11) is indeterminate; 
  // in this case, take σ1= 0
  if (double_almost_equal(lat1, 0e0, 2) && double_almost_equal(a1, DPI/2e0))
    __sigma1 = 0e0;
  // const double sigma1 = std::atan2(sinbeta1, cosalpha * cosbeta1);
  const double sigma1 = __sigma1;
  //printf("\nsigma1=%15.11f", rad2deg(sigma1));
  const double sinsigma1 = std::sin(sigma1);
  const double cossigma1 = std::cos(sigma1);
  const double omega1 = std::atan2(sinalpha0 * sinsigma1, cossigma1);
  //printf("\nomega1=%15.11f", rad2deg(omega1));
  // cmplx_norm =
  // std::sqrt((cosalpha0*cossigma)*(cosalpha0*cossigma)+sinalpha0*sinalpha0);
  // beta = std::atan2(cosalpha0*sinsigma, cmplx_norm);
  // alpha = std::atan2(sinalpha0, cosalpha0*cossigma);

  // Determine sigma2
  const double sec_ecc2 = (a * a - b * b) / (b * b);
  const double cosalpha0 = std::cos(alpha0);
  const double ksqrd = sec_ecc2 * cosalpha0 * cosalpha0;
  //printf("\nk2=%17.14f", ksqrd);
  const double onepk2 = std::sqrt(1e0 + ksqrd);
  const double epsilon = (onepk2 - 1e0) / (onepk2 + 1e0);
  //printf("\nepsilon=%17.14f", epsilon);
  const double ep2 = epsilon * epsilon;
  const double Alpha1 =
      (1e0 + ep2 * (0.25e0 + ep2 * (1e0 / 64e0 + (1e0 / 256e0) * ep2))) /
      (1e0 - epsilon);
  //printf("\nAlpha1=%17.14f", Alpha1);
  const double giota = series_expansion_I(epsilon, sigma1) * Alpha1;
  //printf("\nI1(s1)=%17.14f", giota);
  const double s1 = giota * b;
  //printf("\ns1=%17.6f", s1);
  const double s2 = s1 + s;
  //printf("\ns2=%17.6f", s2);
  const double tau2 = s2 / (b * Alpha1);
  //printf("\ntau2=%17.11f", rad2deg(tau2));
  const double sigma2 = series_expansion_S(epsilon, tau2);
  //printf("\nsigma2=%17.11f", rad2deg(sigma2));

  // Solve triangle NEB
  const double sinsigma2 = std::sin(sigma2);
  const double cossigma2 = std::cos(sigma2);
  const double alpha2 = std::atan2(sinalpha0, cosalpha0 * cossigma2);
  //printf("\nalpha2=%15.11f", rad2deg(alpha2));
  cmplx_norm = std::sqrt((cosalpha0 * cossigma2) * (cosalpha0 * cossigma2) +
                         sinalpha0 * sinalpha0);
  const double beta2 = std::atan2(cosalpha0 * sinsigma2, cmplx_norm);
  //printf("\nbeta2=%15.11f", rad2deg(beta2));
  const double omega2 = std::atan2(sinalpha0 * sinsigma2, cossigma2);
  //printf("\nomega2=%15.11f", rad2deg(omega2));

  // determine λ12
  const double n = third_flattening(f);
  const double A3 =
      1e0 +
      epsilon *
          (-(0.5e0 - 0.5e0 * n) +
           epsilon * (-((1e0 / 4e0) + n * ((1e0 / 8e0) - (3e0 / 8e0) * n)) +
                      epsilon * (-((1e0 / 16e0) +
                                   n * ((3e0 / 16e0) + (1e0 / 16e0) * n)) +
                                 epsilon * (-((3e0 / 64e0) + (1e0 / 32e0) * n) -
                                            (3e0 / 128e0) * epsilon))));
  //printf("\nA3=%17.11f", A3);
  double sigmas[] = {sigma1, sigma2}, giotas[2];
  series_expansion_I3(epsilon, n, sigmas, giotas); // perform both expansions
  const double giota3_s1 = A3*giotas[0];
  //printf("\nI3(s1)=%17.11f", giota3_s1);
  const double giota3_s2 = A3*giotas[1];
  //printf("\nI3(s2)=%17.11f", giota3_s2);
  const double lambda1 = omega1 - f*std::sin(alpha0)*giota3_s1;
  //printf("\nlambda1=%17.11f", rad2deg(lambda1));
  const double lambda2 = omega2 - f*std::sin(alpha0)*giota3_s2;
  //printf("\nlambda2=%17.11f", rad2deg(lambda2));
  const double lambda12 = lambda2 - lambda1;
  //printf("\nlambda12=%17.11f", rad2deg(lambda12));

  // results
  lat2 = std::atan2(std::tan(beta2), (1e0-f));
  //printf("\nlat2=%17.11f", rad2deg(lat2));
  lon2 = ngpt::normalize_angle(lambda2, -ngpt::DPI, ngpt::DPI);
  lon2 = ngpt::normalize_angle(lambda12, -ngpt::DPI, ngpt::DPI) + ngpt::normalize_angle(lon1, -ngpt::DPI, ngpt::DPI);

  return alpha2;
}
