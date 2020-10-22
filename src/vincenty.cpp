#include "vincenty.hpp"
#include "geoconst.hpp"
#include "geodesy.hpp"
#include <cmath>
#include <stdexcept>

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
    const double a1 = 1e0 - f/4e0 -f*f/16e0;
    const double a3 = f/4e0 - f*f/8e0;
    const double a5 = 3e0*f / 16e0;
    const double b1 = 1e0 + f/4e0 +f*f/8e0;
    const double b3 = 1e0 - b1;
    const double b5 = 0e0;
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
  const int MAX_ITERATIONS = 100;
  int iteration = 0;

  const double a = semi_major;
  const double f = flattening;
  const double b = semi_minor;

  double L{lon2 - lon1};
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
  double lambda{L};
  double sinSigma, cosSigma, sigma, sinAlpha, cosSqAlpha, C, cos2SigmaM,
      lambdaP, sinLambda, cosLambda;

  do {
    if (++iteration > MAX_ITERATIONS) {
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

/// @brief Direct Vincenty formula.
///
/// Given an initial point (lat1, lon1), an initial azimuth, a1, and a distance,
/// s, along the geodesic the problem is to find the end point (lat2, lon2)
/// and azimuth, a2.
/// Vincenty’s solution for the distance between points on an ellipsoidal earth
/// model is accurate to within 0.5 mm distance (!), 0.000015'' bearing, on the
/// ellipsoid being used. Calculations based on a spherical earth model, such
/// as the (much simpler) Haversine, are accurate to around 0.3% - which is
/// still good enough for many (most?) purposes, of course (see
/// https://www.movable-type.co.uk/scripts/latlong-vincenty.html).
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
/// @see https://en.wikipedia.org/wiki/Vincenty%27s_formulae
double ngpt::core::direct_vincenty(double lat1, double lon1, double a1,
                                   double s, double semi_major,
                                   double flattening, double semi_minor,
                                   double &lat2, double &lon2,
                                   double convergence_limit) {
  const int MAX_ITERATIONS = 100;
  int iteration = 0;

  double a = semi_major;
  double f = flattening;
  double b = semi_minor;

  const double cosa1{std::cos(a1)};
  const double sina1{std::sin(a1)};
  const double U1{std::atan((1e0 - f) * std::tan(lat1))};
  const double sigma1{std::atan2(std::tan(U1), cosa1)};
  const double sina{std::cos(U1) * sina1};
  const double cosaSq{1e0 - sina * sina};
  const double uSq{cosaSq * (a * a - b * b) / (b * b)};
  const double k1{(std::sqrt(1e0 + uSq) - 1e0) / (std::sqrt(1e0 + uSq) + 1e0)};
  const double A{(1e0 + 0.25e0 * k1 * k1) / (1e0 - k1)};
  const double B{k1 * (1e0 - (3e0 / 8e0) * k1 * k1)};

  double sigma{s / b * A}; // initial guess
  double sigmaP, sigmaM2, cosSigmaM2, deltaSigma, sinSigma;

  do {
    sigmaM2 = 2e0 * sigma1 + sigma;
    cosSigmaM2 = std::cos(sigmaM2);
    sinSigma = std::sin(sigma);
    deltaSigma =
        B * sinSigma *
        (cosSigmaM2 +
         (1e0 / 4e0) * B *
             (std::cos(sigma) * (-1e0 + 2e0 * cosSigmaM2 * cosSigmaM2) -
              (1e0 / 6e0) * B * cosSigmaM2 *
                  (-3e0 + 4e0 * sinSigma * sinSigma) *
                  (-3e0 + 4e0 * cosSigmaM2 * cosSigmaM2)));
    sigmaP = sigma;
    sigma = s / (b * A) + deltaSigma;
  } while (std::abs(sigma - sigmaP) > convergence_limit &&
           ++iteration < MAX_ITERATIONS);
  if (iteration >= MAX_ITERATIONS) {
    throw std::out_of_range(
        "[ERROR] Direct Vincenty cannot converge after 100 iterations!");
  }

  const double sinU1{std::sin(U1)};
  const double cosU1{std::cos(U1)};
  const double cosSigma{std::cos(sigma)};

  // compute latitude
  double nom{sinU1 * cosSigma + cosU1 * sinSigma * cosa1};
  double denom{sina * sina +
               std::pow(sinU1 * sinSigma - cosU1 * cosSigma * cosa1, 2e0)};
  denom = std::sqrt(denom) * (1e0 - f);
  lat2 = std::atan2(nom, denom);

  // compute longtitude
  nom = sinSigma * sina1;
  denom = cosU1 * cosSigma - sinU1 * sinSigma * cosa1;
  double lambda{std::atan2(nom, denom)};
  double C{(f / 16e0) * cosaSq * (4e0 + f * (4e0 - 3e0 * cosaSq))};
  double L{lambda -
           (1e0 - C) * f * sina *
               (sigma +
                C * sinSigma *
                    (cosSigmaM2 +
                     C * cosSigma * (-1e0 + 2e0 * cosSigmaM2 * cosSigmaM2)))};
  // there are some rare case (e.g when the two points run between vertices
  // (i.e. a1 = a2 = 90deg) or end close to vertices, when the following
  // result (lon2) takes values (a little less) than 180degrees (due to roundoff
  // errors. To protect against this, we normalize the longtitude in the range
  // [-180, 180] degrees.
  lon2 = ngpt::normalize_angle(L + lon1, -ngpt::DPI, ngpt::DPI);

  // compute azimouth
  return std::atan2(sina, -sinU1 * sinSigma + cosU1 * cosSigma * cosa1);
}

/// @brief Direct Vincenty formula.
///
/// Given an initial point (lat1, lon1), an initial azimuth, a1, and a distance,
/// s, along the geodesic the problem is to find the end point (lat2, lon2)
/// and azimuth, a2.
/// Vincenty’s solution for the distance between points on an ellipsoidal earth
/// model is accurate to within 0.5 mm distance (!), 0.000015'' bearing, on the
/// ellipsoid being used. Calculations based on a spherical earth model, such
/// as the (much simpler) Haversine, are accurate to around 0.3% - which is
/// still good enough for many (most?) purposes, of course (see
/// https://www.movable-type.co.uk/scripts/latlong-vincenty.html).
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
/// @see https://www.movable-type.co.uk/scripts/latlong-vincenty.html
double ngpt::core::direct_vincenty2(double lat1, double lon1, double a1,
                                    double s, double semi_major,
                                    double flattening, double semi_minor,
                                    double &lat2, double &lon2,
                                    double convergence_limit) {
  const int MAX_ITERATIONS = 100;
  int iteration = 0;

  double a = semi_major;
  double f = flattening;
  double b = semi_minor;

  const double cosa1{std::cos(a1)};
  const double sina1{std::sin(a1)};
  const double tanU1 = (1e0 - f) * std::tan(lat1);
  const double cosU1 = 1e0 / std::sqrt((1e0 + tanU1 * tanU1));
  const double sinU1 = tanU1 * cosU1;
  const double sigma1 = std::atan2(
      tanU1,
      cosa1); // σ1 = angular distance on the sphere from the equator to P1
  const double sina = cosU1 * sina1;
  const double cosSqa = 1e0 - sina * sina;
  const double uSq = cosSqa * (a * a - b * b) / (b * b);
  const double A =
      1e0 +
      uSq / 16384e0 * (4096e0 + uSq * (-768e0 + uSq * (320e0 - 175e0 * uSq)));
  const double B =
      uSq / 1024e0 * (256e0 + uSq * (-128e0 + uSq * (74e0 - 47e0 * uSq)));
  double sigma = s / (b * A), sinSigma = 0e0, cosSigma = 0e0,
         deltaSigma = 0e0;  // σ = angular distance P₁ P₂ on the sphere
  double cos2sigma_m = 0e0; // σₘ = angular distance on the sphere from the
                            // equator to the midpoint of the line
  double sigma_t = 0e0;
  do {
    cos2sigma_m = std::cos(2e0 * sigma1 + sigma);
    sinSigma = std::sin(sigma);
    cosSigma = std::cos(sigma);
    deltaSigma =
        B * sinSigma *
        (cos2sigma_m +
         B / 4e0 *
             (cosSigma * (-1e0 + 2e0 * cos2sigma_m * cos2sigma_m) -
              B / 6e0 * cos2sigma_m * (-3e0 + 4e0 * sinSigma * sinSigma) *
                  (-3e0 + 4e0 * cos2sigma_m * cos2sigma_m)));
    sigma_t = sigma;
    sigma = s / (b * A) + deltaSigma;
  } while (std::abs(sigma - sigma_t) > convergence_limit &&
           ++iteration < MAX_ITERATIONS);

  if (iteration >= MAX_ITERATIONS) {
    throw std::out_of_range(
        "[ERROR] Direct Vincenty cannot converge after 100 iterations!");
  }

  const double x = sinU1 * sinSigma - cosU1 * cosSigma * cosa1;
  lat2 = std::atan2(sinU1 * cosSigma + cosU1 * sinSigma * cosa1,
                    (1e0 - f) * std::sqrt(sina * sina + x * x));
  const double lambda =
      std::atan2(sinSigma * sina1, cosU1 * cosSigma - sinU1 * sinSigma * cosa1);
  const double C = f / 16e0 * cosSqa * (4e0 + f * (4e0 - 3e0 * cosSqa));
  const double L =
      lambda - (1e0 - C) * f * sina *
                   (sigma + C * sinSigma *
                                (cos2sigma_m +
                                 C * cosSigma *
                                     (-1e0 + 2e0 * cos2sigma_m * cos2sigma_m)));
  // there are some rare case (e.g when the two points run between vertices
  // (i.e. a1 = a2 = 90deg) or end close to vertices, when the following
  // result (lon2) takes values (a little less) than 180degrees (due to roundoff
  // errors. To protect against this, we normalize the longtitude in the range
  // [-180, 180] degrees.
  lon2 = ngpt::normalize_angle(L + lon1, -ngpt::DPI, ngpt::DPI);

  const double a2 = std::atan2(sina, -x);
  return std::fmod(a2 + ngpt::D2PI, D2PI);
}
