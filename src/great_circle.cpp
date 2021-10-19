///
/// @file great_circle.cpp
///
/// @todo

#include "great_circle.hpp"
#include "geoconst.hpp"
#include <cassert>

/// @brief Compute the haversine function.
double dso::core::haversine_angle(double angle) noexcept {
  const double sinHalfTheta{std::sin(angle / 2e0)};
  return sinHalfTheta * sinHalfTheta;
}

/// @brief Compute great circle distance between two points on a sphere of
///        radius R
double dso::core::great_circle_distance_cosines(double lat1, double lon1,
                                                 double lat2, double lon2,
                                                 double R) noexcept {
  const double sinf1 = std::sin(lat1);
  const double sinf2 = std::sin(lat2);
  const double cosf1 = std::cos(lat1);
  const double cosf2 = std::cos(lat2);
  const double dlambda = std::abs(lon1 - lon2);
  const double cosdl = std::cos(dlambda);
  // central angle
  const double dsigma = std::acos(sinf1 * sinf2 + cosf1 * cosf2 * cosdl);
  return R * dsigma;
}

/// @brief Compute great circle distance between two points on a sphere of
///        radius R
double dso::core::great_circle_distance_haversine(double lat1, double lon1,
                                                   double lat2, double lon2,
                                                   double R) noexcept {
  const double h{haversine_angle(lat2 - lat1) +
                 std::cos(lat1) * std::cos(lat2) *
                     haversine_angle(lon2 - lon1)};
  return 2e0 * R * std::asin(std::sqrt(h));
}

/// @brief Compute great circle distance between two points on a sphere of
///        radius R
double dso::core::great_circle_distance_vincenty(double lat1, double lon1,
                                                  double lat2, double lon2,
                                                  double R) noexcept {
  const double sinf1 = std::sin(lat1);
  const double sinf2 = std::sin(lat2);
  const double cosf1 = std::cos(lat1);
  const double cosf2 = std::cos(lat2);
  const double dlambda = std::abs(lon1 - lon2);
  const double cosdl = std::cos(dlambda);
  const double sindl = std::sin(dlambda);
  const double c1 = cosf2 * sindl;
  const double c2 = cosf1 * sinf2 - sinf1 * cosf2 * cosdl;
  const double c3 = sinf1 * sinf2 + cosf1 * cosf2 * cosdl;
  // central angle
  const double dsigma = std::atan2(std::sqrt(c1 * c1 + c2 * c2), c3);
  return R * dsigma;
}

/// @brief Compute great circle distance between two points on a sphere of
///        radius R
double dso::core::great_circle_distance_chordl(double lat1, double lon1,
                                                double lat2, double lon2,
                                                double R) noexcept {
  const double sinf1 = std::sin(lat1);
  const double sinf2 = std::sin(lat2);
  const double sinl1 = std::sin(lon1);
  const double sinl2 = std::sin(lon2);
  const double cosf1 = std::cos(lat1);
  const double cosf2 = std::cos(lat2);
  const double cosl1 = std::cos(lon1);
  const double cosl2 = std::cos(lon2);
  const double dx = cosf2 * cosl2 - cosf1 * cosl1;
  const double dy = cosf2 * sinl2 - cosf1 * sinl1;
  const double dz = sinf2 - sinf1;
  const double c = std::sqrt(dx * dx + dy * dy + dz * dz);
  // assert( c>=-1e0 && c<=1e0 );
  // central angle
  const double dsigma = std::asin(c / 2e0) * 2e0;
  return R * dsigma;
}

/// @brief Compute great circle distance between two points on a sphere of
///        radius R
double dso::core::great_circle_distance_pythagoras(double lat1, double lon1,
                                                    double lat2, double lon2,
                                                    double R) noexcept {
  const double fm = (lat1 + lat2) / 2e0;
  const double x = (lon2 - lon1) * std::cos(fm);
  const double y = (lat2 - lat1);
  return std::sqrt(x * x + y * y) * R;
}

/// @brief Compute great circle distance between two points on a sphere of
///        radius R
double dso::core::great_circle_distance_polar(double lat1, double lon1,
                                               double lat2, double lon2,
                                               double R) noexcept {
  const double theta1 = dso::DPI - lat1;
  const double theta2 = dso::DPI - lat2;
  return R * std::sqrt(theta1 * theta1 + theta2 * theta2 -
                       2e0 * theta1 * theta2 * std::cos(lon2 - lon1));
}
