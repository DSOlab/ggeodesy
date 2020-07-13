///
/// @file vincenty.hpp
///
/// @brief Implementation of Vincenty's formulae and other algorithms to
///        compute great circle/godesic quantities.
///
/// Core function/algorithm implementations are nested in the ngpt::core
/// namespace. Users can access the core functions via ngpt functions that
/// either take ellpsoid parameters as compile-time constants (aka template
/// implementations) or runtime parameters (where an Ellipsoid instance is
/// needed).
///
/// @todo
///       - Need to test!
///       - Implement the direct/inverse algorithms from GographicLib, see
///         http://geographiclib.sourceforge.net/html/classGeographicLib_1_1Geodesic.html

#ifndef __NGPT_VINCENTY_HPP__
#define __NGPT_VINCENTY_HPP__

#include "ellipsoid.hpp"
#include <cmath>
#include <stdexcept>

namespace ngpt {

namespace core {

/// @brief Compute the haversine function.
///
/// This function is used within the haversine formula to compute great-circle
/// distances between two points on a sphere from their longitudes and
/// latitudes.
/// It is called the "haversine function" (half a versine) and given an angle
/// θ it is computed as:
/// \f$ hav(\theta ) =  \sin^2 \frac{\theta}{2} = \frac{1- \cos \theta}{2}\f$
///
/// @param[in] angle  The angle to compute the haversine function for (radians)
/// @return The haversine function
///
/// @see https://en.wikipedia.org/wiki/Haversine_formula
double haversine_angle(double angle) noexcept {
  const double sinHalfTheta{std::sin(angle / 2e0)};
  return sinHalfTheta * sinHalfTheta;
}

/// @brief Compute the great circle distance between two points using the
/// haversine formula.
///
/// Given the (ellipsoidal) coordinates of two points (1 and 2) in radians,
/// compute the great circle distance between them, using the haversine formula
/// see https://en.wikipedia.org/wiki/Haversine_formula.
///
/// For this computation we need the radius of the sphere/ellipsoid. Here, we
/// considere the mean radius as: \f$ \frac{2*semi-major+semi-minor}{3} \f$
///
/// @param[in] lat1  Latitude of point 1 (radians)
/// @param[in] lon1  Longtitude  point 1 (radians)
/// @param[in] lat2  Latitude of point 2 (radians)
/// @param[in] lon2  Longtitude  point 2 (radians)
/// @param[in] a     Semi-major axis of reference ellipsoid (meters)
/// @param[in] b     Semi-minor axis of reference ellipsoid (meters)
/// @return          Great circle distance between points 1 and 2 (meters)
///
/// @warning h only approaches 1 for antipodal points (on opposite sides of the
///          sphere). This will cause the formula to fail!
///
/// @note The haversine formula and law of cosines can't be guaranteed correct
///       to better than 0.5%
///
/// @see https://en.wikipedia.org/wiki/Haversine_formula
///
double haversine(double lat1, double lon1, double lat2, double lon2, double a,
                 double b) {
  const double h{haversine_angle(lat2 - lat1) +
                 std::cos(lat1) * std::cos(lat2) *
                     haversine_angle(lon2 - lon1)};
  // according to wikipedia, The International Union of Geodesy and Geophysics
  // (IUGG) defines the mean radius (denoted R_1) to be
  // 2a+b/3
  const double EarthRadius{(2e0 * a + b) / 3e0};
  return 2e0 * EarthRadius * std::asin(std::sqrt(h));
}

/// @brief Compute the inverse Vincenty formula.
///
/// Given the (ellipsoidal) coordinates of two points (1 and 2) in radians,
/// calculate the forward azimouths from 1->2 (i.e. a12) and from 2->1 (i.e.
/// a21) and the ellipsoidal (along the geodesic) distance between the two
/// points. The inverse Vincenty's formula is used for the computation.
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
/// @param[in]  convergence_limit Convergence limit (radians). 10e−12
///                         corresponds to approximately 0.06mm
/// @return      Ellipsoidal (along the geodesic) distance between the two
///              points (meters).
///
/// @throw Will throw an std::out_of_range (implying the input aruments/points)
///        in case more than MAX_ITERATIONS(=100) are needed for the algorithm
///        to converge. Solution to the inverse problem fails to converge or
///        converges slowly for nearly antipodal points.
///
/// @see https://en.wikipedia.org/wiki/Vincenty's_formulae
///
double inverse_vincenty(double lat1, double lon1, double lat2, double lon2,
                        double semi_major, double flattening, double semi_minor,
                        double &a12, double &a21,
                        double convergence_limit = 1e-12) {
  const int MAX_ITERATIONS = 100;
  int iteration = 0;

  const double a = semi_major;
  const double f = flattening;
  const double b = semi_minor;

  const double U1{std::atan((1e0 - f) * std::tan(lat1))};
  const double U2{std::atan((1e0 - f) * std::tan(lat2))};
  const double L{lon2 - lon1};
  const double sinU1{std::sin(U1)};
  const double sinU2{std::sin(U2)};
  const double cosU1{std::cos(U1)};
  const double cosU2{std::cos(U2)};
  double lambda{L};
  double sinSigma, cosSigma, sigma, sinAlpha, cosSqAlpha, C, cos2SigmaM,
      lambdaP, sinLambda, cosLambda;

  do {
    if (++iteration > MAX_ITERATIONS) {
      throw std::out_of_range(
          "[ERROR] Inverse Vincenty cannot converge after 100 iterations!");
    }
    sinLambda = std::sin(lambda);
    cosLambda = std::cos(lambda);
    sinSigma = std::sqrt((cosU2 * sinLambda) * (cosU2 * sinLambda) +
                         (cosU1 * sinU2 - sinU1 * cosU2 * cosLambda) *
                             (cosU1 * sinU2 - sinU1 * cosU2 * cosLambda));
    cosSigma = sinU1 * sinU2 + cosU1 * cosU2 * cosLambda;
    sigma = atan2(sinSigma, cosSigma);
    sinAlpha = cosU1 * cosU2 * sinLambda / sinSigma;
    cosSqAlpha = 1e0 - sinAlpha * sinAlpha;
    cos2SigmaM = cosSigma - 2e0 * sinU1 * sinU2 / cosSqAlpha;
    C = (f / 16e0) * cosSqAlpha * (4e0 + f * (4e0 - 3e0 * cosSqAlpha));
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
/// @param[in]  convergence_limit Convergence limit (radians). 10e−12
///                         corresponds to approximately 0.06mm
/// @return     Azimouth from point 2 to point 1 (radians)
///
/// @see https://en.wikipedia.org/wiki/Vincenty%27s_formulae
double direct_vincenty(double lat1, double lon1, double a1, double s,
                       double semi_major, double flattening, double semi_minor,
                       double &lat2, double &lon2,
                       double convergence_limit = 1e-12) {
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
    if (++iteration > MAX_ITERATIONS) {
      throw std::out_of_range(
          "[ERROR] Direct Vincenty cannot converge after 100 iterations!");
    }
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
  } while (std::abs(sigma - sigmaP) > convergence_limit);

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
  lon2 = L + lon1;

  // compute azimouth
  return std::atan2(sina, -sinU1 * sinSigma + cosU1 * cosSigma * cosa1) + DPI;
}

} // namespace core

/// @brief Compute the great circle distance between two points using the
/// haversine formula.
///
/// @tparam    E     ngpt::ellipsoid enum
/// @param[in] lat1  Latitude of point 1 (radians)
/// @param[in] lon1  Longtitude  point 1 (radians)
/// @param[in] lat2  Latitude of point 2 (radians)
/// @param[in] lon2  Longtitude  point 2 (radians)
/// @return          Great circle distance between points 1 and 2 (meters)
///
/// @see ngpt::core::haversine
template <ellipsoid E = ellipsoid::wgs84>
double haversine(double lat1, double lon1, double lat2, double lon2) {
  const double a = ellipsoid_traits<E>::a;
  const double b = semi_minor<E>();
  return core::haversine(lat1, lon1, lat2, lon2, a, b);
}

/// @brief Compute the great circle distance between two points using the
/// haversine formula.
///
/// @param[in] e     Reference Ellipsoid
/// @param[in] lat1  Latitude of point 1 (radians)
/// @param[in] lon1  Longtitude  point 1 (radians)
/// @param[in] lat2  Latitude of point 2 (radians)
/// @param[in] lon2  Longtitude  point 2 (radians)
/// @return          Great circle distance between points 1 and 2 (meters)
///
/// @see ngpt::core::haversine
double haversine(double lat1, double lon1, double lat2, double lon2,
                 const Ellipsoid &e) {
  const double a = e.semi_major();
  const double b = e.semi_minor();
  return core::haversine(lat1, lon1, lat2, lon2, a, b);
}

/// @brief Compute the inverse Vincenty formula.
///
/// @tparam     E           ngpt::ellipsoid enum
/// @param[in]  lat1        Latitude of point 1 (radians)
/// @param[in]  lon1        Longtitude of point 1 (radians)
/// @param[in]  lat2        Latitude of point 2 (radians)
/// @param[in]  lon2        Longtitude of point 2 (radians)
/// @param[out] a12         Azimouth from point 1 to point 2 (radians)
/// @param[out] a21         Azimouth from point 2 to point 1 (radians)
/// @param[in]  convergence_limit Convergence limit (radians). 10e−12
///                         corresponds to approximately 0.06mm
/// @return      Ellipsoidal (along the geodesic) distance between the two
///              points (meters).
///
/// @see ngpt::core::inverse_vincenty
template <ellipsoid E = ellipsoid::wgs84>
double inverse_vincenty(double lat1, double lon1, double lat2, double lon2,
                        double &a12, double &a21,
                        double convergence_limit = 1e-12) {
  const double a = ellipsoid_traits<E>::a;
  const double f = ellipsoid_traits<E>::f;
  const double b = semi_minor<E>();
  return core::inverse_vincenty(lat1, lon1, lat2, lon2, a, f, b, a12, a21,
                                convergence_limit);
}

/// @brief Compute the inverse Vincenty formula.
///
/// @param[in]  e           Reference Ellipsoid
/// @param[in]  lat1        Latitude of point 1 (radians)
/// @param[in]  lon1        Longtitude of point 1 (radians)
/// @param[in]  lat2        Latitude of point 2 (radians)
/// @param[in]  lon2        Longtitude of point 2 (radians)
/// @param[out] a12         Azimouth from point 1 to point 2 (radians)
/// @param[out] a21         Azimouth from point 2 to point 1 (radians)
/// @param[in]  convergence_limit Convergence limit (radians). 10e−12
///                         corresponds to approximately 0.06mm
/// @return      Ellipsoidal (along the geodesic) distance between the two
///              points (meters).
///
/// @see ngpt::core::inverse_vincenty
double inverse_vincenty(double lat1, double lon1, double lat2, double lon2,
                        double &a12, double &a21, const Ellipsoid &e,
                        double convergence_limit = 1e-12) {
  const double a = e.semi_major();
  const double f = e.flattening();
  const double b = e.semi_minor();
  return core::inverse_vincenty(lat1, lon1, lat2, lon2, a, f, b, a12, a21,
                                convergence_limit);
}

/// @brief Direct Vincenty formula.
///
/// @tparam     E           ngpt::ellpsoid enum
/// @param[in]  lat1        Latitude of point 1 (radians)
/// @param[in]  lon1        Longtitude of point 1 (radians)
/// @param[in]  a1          Azimouth from point 1 to point 2 (radians)
/// @param[in]  s           Ellipsoid distance (along the geodesic) from point
///                         1 to point 2 (meters)
/// @param[out] lat2        Latitude of point 2 (radians)
/// @param[out] lon2        Longtitude of point 2 (radians)
/// @param[in]  convergence_limit Convergence limit (radians). 10e−12
///                         corresponds to approximately 0.06mm
/// @return     Azimouth from point 2 to point 1 (radians)
///
/// @see ngpt::core::direct_vincenty
template <ellipsoid E = ellipsoid::wgs84>
double direct_vincenty(double lat1, double lon1, double a1, double s,
                       double &lat2, double &lon2,
                       double convergence_limit = 1e-12) {
  const double a = ellipsoid_traits<E>::a;
  const double f = ellipsoid_traits<E>::f;
  const double b = semi_minor<E>();
  return core::direct_vincenty(lat1, lon1, a1, s, a, f, b, lat2, lon2,
                               convergence_limit);
}
/// @brief Direct Vincenty formula.
///
/// @param[in]  e           Reference Ellipsoid
/// @param[in]  lat1        Latitude of point 1 (radians)
/// @param[in]  lon1        Longtitude of point 1 (radians)
/// @param[in]  a1          Azimouth from point 1 to point 2 (radians)
/// @param[in]  s           Ellipsoid distance (along the geodesic) from point
///                         1 to point 2 (meters)
/// @param[out] lat2        Latitude of point 2 (radians)
/// @param[out] lon2        Longtitude of point 2 (radians)
/// @param[in]  convergence_limit Convergence limit (radians). 10e−12
///                         corresponds to approximately 0.06mm
/// @return     Azimouth from point 2 to point 1 (radians)
///
/// @see ngpt::core::direct_vincenty
double direct_vincenty(double lat1, double lon1, double a1, double s,
                       double &lat2, double &lon2, const Ellipsoid &e,
                       double convergence_limit = 1e-12) {
  const double a = e.semi_major();
  const double f = e.flattening();
  const double b = e.semi_minor();
  return core::direct_vincenty(lat1, lon1, a1, s, a, f, b, lat2, lon2,
                               convergence_limit);
}

} // end namespace ngpt
#endif
