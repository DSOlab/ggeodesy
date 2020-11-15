#include "vincenty.hpp"
#include "geoconst.hpp"
#include <cmath>
#include <stdexcept>
#ifdef DEBUG
#include <fenv.h>
#endif

using ngpt::D2PI;
using ngpt::DPI;

/// @brief Compute the inverse Vincenty formula.
///
/// Given the (ellipsoidal) coordinates of two points (1 and 2) in radians,
/// calculate the forward azimouths from 1->2 (i.e. a12) and from 2->1 (i.e.
/// a21) and the ellipsoidal (along the geodesic) distance between the two
/// points. The inverse Vincenty's formula is used for the computation.
/// Vincenty’s solution for the distance between points on an ellipsoidal earth
/// model is accurate to within 0.5 mm distance (!), 0.000015'' bearing, on the
/// ellipsoid being used. Calculations based on a spherical earth model, such
/// as the (much simpler) Haversine, are accurate to around 0.3% – which is
/// still good enough for many (most?) purposes, of course (see
/// https://www.movable-type.co.uk/scripts/latlong-vincenty.html).
///
/// @param[in]  lat1        Latitude of point 1 (radians)
/// @param[in]  lon1        Longtitude of point 1 (radians)
/// @param[in]  lat2        Latitude of point 2 (radians)
/// @param[in]  lon2        Longtitude of point 2 (radians)
/// @param[in]  semi_major  Semi-major axis of reference ellipsoid (meters)
/// @param[in]  semi_minor  Semi-minor axis of reference ellipsoid (meters)
/// @param[in]  flattening  Flattening of reference ellipsoid
/// @param[out] a12         Azimouth from point 1 to point 2 (radians)
/// @param[out] a21         Azimouth from point 2 to point 1 (radians)
/// @param[in]  convergence_limit Convergence limit (radians). 10e-12
///                         corresponds to approximately 0.06mm
/// @return      Ellipsoidal (along the geodesic) distance between the two
///              points (meters).
///
/// @throw Will throw an std::out_of_range (implying the input aruments/points)
///        in case more than MAX_ITERATIONS(=100) are needed for the algorithm
///        to converge. Solution to the inverse problem fails to converge or
///        converges slowly for nearly antipodal points.
///
/// @see [1] https://en.wikipedia.org/wiki/Vincenty's_formulae
/// @see [2] https://futureboy.us/fsp/colorize.fsp?fileName=navigation.frink
/// @see [3] https://github.com/boostorg/geometry/blob/develop/include/boost/geometry/formulas/vincenty_inverse.hpp
/// @see [4] https://geographiclib.sourceforge.io/geodesic-papers/vincenty75b.pdf
double ngpt::core::inverse_vincenty_nearly_antipodal(
    double lat1, double lon1, double lat2, double lon2, double semi_major,
    double flattening, double semi_minor, double &a12, double &a21,
    double convergence_limit) {
  const int MAX_ITERATIONS = 1000;
  int iteration = 0;

  const double a = semi_major;
  const double f = flattening;
  const double b = semi_minor;

  double L{lon2 - lon1};
  if (L < -DPI)
    L += D2PI;
  if (L > DPI)
    L -= D2PI;
  const double Ld = (L > 0e0) ? (DPI - L) : (-DPI - L);
  //  For the trigonometrics of U1 and U2, boost uses the following:
  //+ see Reference [3]
  const double one_min_f = 1e0 - f;
  const double tanU1 = one_min_f * std::tan(lat1);
  const double tanU2 = one_min_f * std::tan(lat2);
  // calculate sinU and cosU using trigonometric identities
  const double temp_den_U1 = std::sqrt(1e0 + std::sqrt(tanU1));
  const double temp_den_U2 = std::sqrt(1e0 + std::sqrt(tanU2));
  // cos = 1 / sqrt(1 + tan^2)
  const double cosU1 = 1e0 / temp_den_U1;
  const double cosU2 = 1e0 / temp_den_U2;
  // sin = tan / sqrt(1 + tan^2)
  // sin = tan * cos
  const double sinU1 = tanU1 * cosU1;
  const double sinU2 = tanU2 * cosU2;
  /*
  const double U1{std::atan((1e0 - f) * std::tan(lat1))};
  const double U2{std::atan((1e0 - f) * std::tan(lat2))};
  const double sinU1{std::sin(U1)};
  const double sinU2{std::sin(U2)};
  const double cosU1{std::cos(U1)};
  const double cosU2{std::cos(U2)};
  */

  double lambdad = 0e0;
  double sinSigma, cosSigma, sigma, sinAlpha(0e0), sinAlpha_prev, cosSqAlpha, C,
      cos2SigmaM, sinLambda, cosLambda, D, sinSqSigma;
  cosSqAlpha = 0.5e0;
  cos2SigmaM = 0e0;
  sigma = DPI - std::abs(std::atan(sinU1 / cosU1) + std::atan(sinU2 / cosU2));
  // special case when φ2 = -φ1
  if (lat2 == -lat1) {
    sigma = DPI;
    lambdad = 0e0;
    //const double a1 = 1e0 - f/4e0 -f*f/16e0;
    //const double a3 = f/4e0 - f*f/8e0;
    //const double a5 = 3e0*f / 16e0;
    const double b1 = 1e0 + f/4e0 +f*f/8e0;
    const double b3 = 1e0 - b1;
    const double b5 = 0e0;
    const double Q = Ld / (f*DPI);
    const double sina = b1*Q+b3*Q*Q*Q + b5*Q*Q*Q*Q*Q;
    a12 = std::asin(sina);
    a12 = std::fmod(a12 + D2PI, D2PI);
    const double epsilon = (a * a - b * b) / (b * b);
    const double E = std::sqrt(1e0 + epsilon * cosSqAlpha);
    const double F = (E - 1e0) / (E + 1e0);
    const double A = (1e0 + 0.25e0 * F * F) / (1e0 - F);
    return b*A*DPI;
  }
  do {
    if (++iteration > MAX_ITERATIONS) {
      auto itstr = std::to_string(iteration);
      throw std::out_of_range(
          "[ERROR] ngpt::core::inverse_vincenty_nearly_antipodal cannot "
          "converge after " +
          itstr + " iterations!");
    }
    sinAlpha_prev = sinAlpha;
    cosSigma = std::cos(sigma);
    sinSigma = std::sin(sigma);
    sinLambda = std::sin(lambdad);
    cosLambda = std::cos(lambdad);
    C = (f / 16e0) * cosSqAlpha * (4e0 + f * (4e0 - 3e0 * cosSqAlpha));
    cos2SigmaM = cosSigma - 2e0 * sinU1 * sinU2 / cosSqAlpha;
    D = (1e0 - C) * f *
        (sigma + C * sinSigma *
                     (cos2SigmaM +
                      C * cosSigma * (-1e0 + 2e0 * cos2SigmaM * cos2SigmaM)));
    sinAlpha = (Ld - lambdad) / D;
    sinLambda = sinAlpha * sinSigma / (cosU1 * cosU2);
    sinSqSigma = (cosU2 * sinLambda) * (cosU2 * sinLambda) +
                 (cosU1 * sinU2 + sinU1 * cosU2 * cosLambda) *
                     (cosU1 * sinU2 + sinU1 * cosU2 * cosLambda);
    cosSqAlpha = 1e0 - sinAlpha * sinAlpha;
    lambdad = std::asin(sinLambda);
    double cosSqSigma = 1e0 - sinSqSigma;
    sigma = std::atan2(std::sqrt(sinSqSigma), std::sqrt(cosSqSigma));
    cosSqAlpha = 1e0 - sinAlpha * sinAlpha;
  } while (std::abs(sinAlpha_prev - sinAlpha) > convergence_limit);

  double SinAlpha1 = sinAlpha / cosU1;
  double CosAlpha1 = std::sqrt(1e0 - SinAlpha1 * SinAlpha1);
  a12 = std::atan2(SinAlpha1, CosAlpha1);
  a12 = std::fmod(a12 + D2PI, D2PI);

  if (cosU1 * sinU2 + sinU1 * cosU2 * std::cos(lambdad))
    CosAlpha1 = std::copysign(CosAlpha1, -1e0);
  else
    CosAlpha1 = std::copysign(CosAlpha1, 1e0);
  a21 = std::atan2(sinAlpha, -sinU1 * sinSigma + cosU1 * cosSigma * CosAlpha1);
  a21 = std::fmod(a21 + D2PI, D2PI);

  const double epsilon = (a * a - b * b) / (b * b);
  const double E = std::sqrt(1e0 + epsilon * cosSqAlpha);
  const double F = (E - 1e0) / (E + 1e0);
  const double A = (1e0 + 0.25e0 * F * F) / (1e0 - F);
  const double B = F * (1e0 - (3e0 / 8e0) * F * F);
  const double deltaSigma =
      B * sinSigma *
      (cos2SigmaM +
       0.25e0 * B *
           (cosSigma * (-1e0 + 2e0 * cos2SigmaM * cos2SigmaM) -
            (1e0 / 6e0) * B * cos2SigmaM * (-3e0 + 4e0 * sinSqSigma) *
                (-3e0 + 4e0 * cos2SigmaM * cos2SigmaM)));
  return b * A * (sigma - deltaSigma);
}

/// @brief Compute the inverse Vincenty formula.
///
/// Given the (ellipsoidal) coordinates of two points (1 and 2) in radians,
/// calculate the forward azimouths from 1->2 (i.e. a12) and from 2->1 (i.e.
/// a21) and the ellipsoidal (along the geodesic) distance between the two
/// points. The inverse Vincenty's formula is used for the computation.
/// Vincenty’s solution for the distance between points on an ellipsoidal earth
/// model is accurate to within 0.5 mm distance (!), 0.000015'' bearing, on the
/// ellipsoid being used. Calculations based on a spherical earth model, such
/// as the (much simpler) Haversine, are accurate to around 0.3% – which is
/// still good enough for many (most?) purposes, of course (see
/// https://www.movable-type.co.uk/scripts/latlong-vincenty.html).
///
/// @param[in]  lat1        Latitude of point 1 (radians)
/// @param[in]  lon1        Longtitude of point 1 (radians)
/// @param[in]  lat2        Latitude of point 2 (radians)
/// @param[in]  lon2        Longtitude of point 2 (radians)
/// @param[in]  semi_major  Semi-major axis of reference ellipsoid (meters)
/// @param[in]  semi_minor  Semi-minor axis of reference ellipsoid (meters)
/// @param[in]  flattening  Flattening of reference ellipsoid
/// @param[out] a12         Azimouth from point 1 to point 2 (radians)
/// @param[out] a21         Azimouth from point 2 to point 1 (radians)
/// @param[in]  convergence_limit Convergence limit (radians). 10e-12
///                         corresponds to approximately 0.06mm
/// @return      Ellipsoidal (along the geodesic) distance between the two
///              points (meters).
///
/// @throw Will throw an std::out_of_range (implying the input aruments/points)
///        in case more than MAX_ITERATIONS(=100) are needed for the algorithm
///        to converge. Solution to the inverse problem fails to converge or
///        converges slowly for nearly antipodal points.
///
/// @see [1] https://en.wikipedia.org/wiki/Vincenty's_formulae
/// @see [2] https://futureboy.us/fsp/colorize.fsp?fileName=navigation.frink
/// @see [3]
/// https://github.com/boostorg/geometry/blob/develop/include/boost/geometry/formulas/vincenty_inverse.hpp
double ngpt::core::inverse_vincenty(double lat1, double lon1, double lat2,
                                    double lon2, double semi_major,
                                    double flattening, double semi_minor,
                                    double &a12, double &a21,
                                    double convergence_limit) {
  const int MAX_ITERATIONS = 1000;
  int iteration = 0;

  const double a = semi_major;
  const double f = flattening;
  const double b = semi_minor;

  double L{lon2 - lon1};
  double lambda{L};
  if (L < -DPI)
    L += D2PI;
  if (L > DPI)
    L -= D2PI;
  //  For the trigonometrics of U1 and U2, boost uses the following:
  //+ see Reference [3]
  const double one_min_f = 1e0 - f;
  const double tanU1 = one_min_f * std::tan(lat1);
  const double tanU2 = one_min_f * std::tan(lat2);
  // calculate sin U and cos U using trigonometric identities
  const double temp_den_U1 = std::sqrt(1e0 + tanU1*tanU1);
  const double temp_den_U2 = std::sqrt(1e0 + tanU2*tanU2);
  // cos = 1 / sqrt(1 + tan^2)
  const double cosU1 = 1e0 / temp_den_U1;
  const double cosU2 = 1e0 / temp_den_U2;
  // sin = tan / sqrt(1 + tan^2)
  // sin = tan * cos
  const double sinU1 = tanU1 * cosU1;
  const double sinU2 = tanU2 * cosU2;
  /*
  const double U1{std::atan((1e0 - f) * std::tan(lat1))};
  const double U2{std::atan((1e0 - f) * std::tan(lat2))};
  const double sinU1{std::sin(U1)};
  const double sinU2{std::sin(U2)};
  const double cosU1{std::cos(U1)};
  const double cosU2{std::cos(U2)};
  */
  double sinSigma, cosSigma, sigma, sinAlpha, cosSqAlpha, C, cos2SigmaM,
      lambdaP, sinLambda, cosLambda;

  do {
    if (++iteration > MAX_ITERATIONS || std::abs(lambda)> DPI) {
      /*try {
        return inverse_vincenty_nearly_antipodal(
            lat1, lon1, lat2, lon2, semi_major, flattening, semi_minor, a12,
            a21, convergence_limit);
      } catch (std::out_of_range &e) {
        printf("\n[ERROR] Failed to converge for nearly antipodal!");
      }*/
      auto itstr = std::to_string(iteration);
      throw std::out_of_range(
          "[ERROR] ngpt::core::inverse_vincenty cannot "
          "converge after " +
          itstr + " iterations!");
    }
    sinLambda = std::sin(lambda);
    cosLambda = std::cos(lambda);
    sinSigma = std::sqrt((cosU2 * sinLambda) * (cosU2 * sinLambda) +
                         (cosU1 * sinU2 - sinU1 * cosU2 * cosLambda) *
                             (cosU1 * sinU2 - sinU1 * cosU2 * cosLambda));
    cosSigma = sinU1 * sinU2 + cosU1 * cosU2 * cosLambda;
    sinAlpha = cosU1 * cosU2 * sinLambda / sinSigma;
    cosSqAlpha = 1e0 - sinAlpha * sinAlpha;
    // check for equatorial points -- see reference [2]
    if (cosSqAlpha == 0e0)
      cos2SigmaM = 0e0;
    else
      cos2SigmaM = cosSigma - 2e0 * sinU1 * sinU2 / cosSqAlpha;
    C = (f / 16e0) * cosSqAlpha * (4e0 + f * (4e0 - 3e0 * cosSqAlpha));
    sigma = atan2(sinSigma, cosSigma);
    lambdaP = lambda;
    lambda = L + (1e0 - C) * f * sinAlpha *
                     (sigma + C * sinSigma *
                                  (cos2SigmaM +
                                   C * cosSigma *
                                       (-1e0 + 2e0 * cos2SigmaM * cos2SigmaM)));
/*
#ifdef DEBUG
  bool isNan = false;
  if (sinLambda != sinLambda) {
    printf("\n[ERROR] sinLambda is Nan");
    isNan = true;
  }
  if (cosLambda != cosLambda) {
    printf("\n[ERROR] cosLambda is Nan");
    isNan = true;
  }
  if (sinSigma != sinSigma) {
    printf("\n[ERROR] sinSigma is Nan");
    isNan = true;
  }
  if (cosSigma != cosSigma) {
    printf("\n[ERROR] cosSigma is Nan");
    isNan = true;
  }
  if (sinAlpha != sinAlpha) {
    printf("\n[ERROR] sinAlpha is Nan");
    isNan = true;
  }
  if (cosSqAlpha != cosSqAlpha) {
    printf("\n[ERROR] cosSqAlpha is Nan");
    isNan = true;
  }
  if (lambda != lambda) {
    printf("\n[ERROR] Lambda is Nan");
    isNan = true;
  }
  if (isNan) {
      auto itstr = std::to_string(iteration);
      throw std::out_of_range(
          "[ERROR] ngpt::core::inverse_vincenty cannot "
          "converge after " +
          itstr + " iterations!");
  }
#endif
  */
  } while (std::abs(lambda - lambdaP) > convergence_limit);

  const double uSq{cosSqAlpha * ((a * a - b * b) / (b * b))};
  const double k1{(std::sqrt(1e0 + uSq) - 1e0) / (std::sqrt(1e0 + uSq) + 1e0)};
  const double A{(1e0 + 0.25e0 * k1 * k1) / (1e0 - k1)};
  const double B{k1 * (1e0 - (3e0 / 8e0) * k1 * k1)};
  const double deltaSigma{
      B * sinSigma *
      (cos2SigmaM +
       B / 4e0 *
           (cosSigma * (-1e0 + 2e0 * cos2SigmaM * cos2SigmaM) -
            B / 6e0 * cos2SigmaM * (-3e0 + 4e0 * sinSigma * sinSigma) *
                (-3e0 + 4e0 * cos2SigmaM * cos2SigmaM)))};
  double distance{b * A * (sigma - deltaSigma)};

  // forward azimouth
  a12 =
      std::atan2(cosU2 * sinLambda, cosU1 * sinU2 - sinU1 * cosU2 * cosLambda);
  // normalize
  a12 = std::fmod(a12 + D2PI, D2PI);

  // backward azimouth
  a21 =
      std::atan2(cosU1 * sinLambda, -sinU1 * cosU2 + cosU1 * sinU2 * cosLambda);
  // normalize
  a21 = std::fmod(a21 + D2PI, D2PI);

  return distance;
}
