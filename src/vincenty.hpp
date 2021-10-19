///
/// @file vincenty.hpp
///
/// @brief Implementation of Vincenty's formulae and other algorithms to
///        compute great circle/godesic quantities.
/// Vincenty's formulae are two related iterative methods used in geodesy to
/// calculate the distance between two points on the surface of a spheroid,
/// developed by Thaddeus Vincenty (1975a). They are based on the assumption
/// that the figure of the Earth is an oblate spheroid, and hence are more
/// accurate than methods that assume a spherical Earth, such as great-circle
/// distance.
/// The first (direct) method computes the location of a point that is a given
/// distance and azimuth (direction) from another point. The second (inverse)
/// method computes the geographical distance and azimuth between two given
/// points. They have been widely used in geodesy because they are accurate to
/// within 0.5 mm (0.020 in) on the Earth ellipsoid.
///
/// Core function/algorithm implementations are nested in the dso::core
/// namespace. Users can access the core functions via ngpt functions that
/// either take ellpsoid parameters as compile-time constants (aka template
/// implementations) or runtime parameters (where an Ellipsoid instance is
/// needed).
///
/// @see
/// https://adl1995.github.io/inaccuracy-in-boost-geometry-geodesic-algorithms-for-nearly-antipodal-points.html
///
/// @todo
///       - Need to test!
///       - Implement the direct/inverse algorithms from GographicLib, see
///         http://geographiclib.sourceforge.net/html/classGeographicLib_1_1Geodesic.html

#ifndef __NGPT_VINCENTY_HPP__
#define __NGPT_VINCENTY_HPP__

#include "ellipsoid.hpp"

namespace dso {

namespace core {

/// @brief Compute the inverse Vincenty formula for nearly antipodal points
double inverse_vincenty_nearly_antipodal(double lat1, double lon1, double lat2,
                                         double lon2, double semi_major,
                                         double flattening, double semi_minor,
                                         double &a12, double &a21,
                                         double convergence_limit = 1e-12);
/// @brief Compute the inverse Vincenty formula.
double inverse_vincenty(double lat1, double lon1, double lat2, double lon2,
                        double semi_major, double flattening, double semi_minor,
                        double &a12, double &a21,
                        double convergence_limit = 1e-12);
double inverse_karney(double lat1, double lon1, double lat2, double lon2,
                        double semi_major, double flattening, double semi_minor,
                        double &a12, double &a21);

/// @brief Direct Vincenty formula.
double direct_vincenty(double lat1, double lon1, double a1, double s,
                       double semi_major, double flattening, double semi_minor,
                       double &lat2, double &lon2,
                       double convergence_limit = 1e-12);
double direct_karney(double lat1, double lon1, double a1, double s,
                     double semi_major, double flattening, double semi_minor,
                     double &lat2, double &lon2);
double direct_vincenty2(double lat1, double lon1, double a1, double s,
                        double semi_major, double flattening, double semi_minor,
                        double &lat2, double &lon2,
                        double convergence_limit = 1e-12);
} // namespace core

/// @brief Compute the inverse Vincenty formula.
///
/// @tparam     E           dso::ellipsoid enum
/// @param[in]  lat1        Latitude of point 1 (radians)
/// @param[in]  lon1        Longtitude of point 1 (radians)
/// @param[in]  lat2        Latitude of point 2 (radians)
/// @param[in]  lon2        Longtitude of point 2 (radians)
/// @param[out] a12         Azimouth from point 1 to point 2 (radians)
/// @param[out] a21         Azimouth from point 2 to point 1 (radians)
/// @param[in]  convergence_limit Convergence limit (radians) \f$10^{-12}\f$
///                         corresponds to approximately 0.06mm
/// @return      Ellipsoidal (along the geodesic) distance between the two
///              points (meters).
///
/// @see dso::core::inverse_vincenty
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
/// @param[in]  convergence_limit Convergence limit (radians) \f$10^{-12}\f$
///                         corresponds to approximately 0.06mm
/// @return      Ellipsoidal (along the geodesic) distance between the two
///              points (meters).
///
/// @see dso::core::inverse_vincenty
/*double inverse_vincenty(double lat1, double lon1, double lat2, double lon2,
                        double &a12, double &a21, const Ellipsoid &e,
                        double convergence_limit = 1e-12) {
  const double a = e.semi_major();
  const double f = e.flattening();
  const double b = e.semi_minor();
  return core::inverse_vincenty(lat1, lon1, lat2, lon2, a, f, b, a12, a21,
                                convergence_limit);
}*/

/// @brief Direct Vincenty formula.
///
/// @tparam     E           dso::ellpsoid enum
/// @param[in]  lat1        Latitude of point 1 (radians)
/// @param[in]  lon1        Longtitude of point 1 (radians)
/// @param[in]  a1          Azimouth from point 1 to point 2 (radians)
/// @param[in]  s           Ellipsoid distance (along the geodesic) from point
///                         1 to point 2 (meters)
/// @param[out] lat2        Latitude of point 2 (radians)
/// @param[out] lon2        Longtitude of point 2 (radians)
/// @param[in]  convergence_limit Convergence limit (radians) \f$10^{-12}\f$
///                         corresponds to approximately 0.06mm
/// @return     Azimouth from point 2 to point 1 (radians)
///
/// @see dso::core::direct_vincenty
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
/// @param[in]  convergence_limit Convergence limit (radians) \f$10^{-12}\f$
///                         corresponds to approximately 0.06mm
/// @return     Azimouth from point 2 to point 1 (radians)
///
/// @see dso::core::direct_vincenty
/*double direct_vincenty(double lat1, double lon1, double a1, double s,
                       double &lat2, double &lon2, const Ellipsoid &e,
                       double convergence_limit = 1e-12) {
  const double a = e.semi_major();
  const double f = e.flattening();
  const double b = e.semi_minor();
  return core::direct_vincenty(lat1, lon1, a1, s, a, f, b, lat2, lon2,
                                convergence_limit);
}*/

} // dso

#endif
