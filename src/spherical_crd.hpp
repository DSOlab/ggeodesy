/// @file spherical_crd.hpp
///
/// @brief Definition and functions for spherical coordinates
///
/// Any spherical coordinate triplet ( r, θ, φ) specifies a single point of
/// three-dimensional space. On the other hand, every point has infinitely many
/// equivalent spherical coordinates. One can add or subtract any number of
/// full turns to either angular measure without changing the angles themselves,
/// and therefore without changing the point. It is also convenient, in many
/// contexts, to allow negative radial distances, with the convention that
/// (−r, θ, φ) is equivalent to (r, θ, φ) for any r, θ, and φ. Moreover,
/// (r, −θ, φ) is equivalent to (r, θ, φ+180).
///
/// If it is necessary to define a unique set of spherical coordinates for
/// each point, one must restrict their ranges. A common choice is
/// * r >= 0,
/// * 0 <= θ <= 180
/// * -180 <= φ <= 180 (sometimes denoted as λ)

#ifndef __SPHERICAL_COORDINATES_NGPT_HPP__
#define __SPHERICAL_COORDINATES_NGPT_HPP__

#include <cmath>

namespace dso {

void spherical2cart(double r, double phi, double theta, double &x, double &y,
                    double &z) noexcept {
  const double sint = std::sin(theta);
  x = r * std::cos(phi) * sint;
  y = r * std::sin(phi) * sint;
  z = r * std::cos(theta);
  return;
}

void cart2spherical(double x, double y, double z, double &r, double &phi,
                    double &theta) noexcept {
  const double r_xy = std::sqrt(x * x + y * y);
  r = std::sqrt(x * x + y * y + z * z);
  // theta = std::atan2(r_xy, z);
  theta = std::acos(z / r);
  phi = std::atan2(y, x);
  return;
}

} // dso
#endif
