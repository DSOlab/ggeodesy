///
/// \file ell2car.hpp
///
/// \brief Transform ellipsoidal to cartesian coordinates.
///

#ifndef __ELLIPSOIDAL_TO_CARTESIAN__
#define __ELLIPSOIDAL_TO_CARTESIAN__

#include <cmath>
#include "ellipsoid.hpp"

namespace ngpt
{

/// @brief Ellipsoidal to cartesian coordinates.
///
/// Transform (geocentric) cartesian coordinates (on the ellipsoid) to
/// ellipsoidal coordinates. Units are meters and radians.
///
/// @tparam     E       The reference ellipsoid (i.e. one of ngpt::ellipsoid).
/// @param[in]   phi    Ellipsoidal latitude, radians.
/// @param[in]   lambda Ellipsoidal longtitude, radians.
/// @param[in]   h      Ellipsoidal height, meters.
/// @param[out]  x      Cartesian x-component, meters.
/// @param[out]  y      Cartesian y-component, meters.
/// @param[out]  z      Cartesian z-component, meters.
/// @throw              Does not throw.
///
template<ellipsoid E>
    void ell2car(double phi, double lambda, double h, double& x, double& y,
        double& z) noexcept
{
    // Eccentricity squared.
    constexpr double e2 { ngpt::eccentricity_squared<E>() };
    
    // Radius of curvature in the prime vertical.
    double N { ngpt::N<E>(phi) };

    // Trigonometric numbers.
    double sinf { std::sin(phi) };
    double cosf { std::cos(phi) };
    double sinl { std::sin(lambda) };
    double cosl { std::cos(lambda) };

    // Compute geocentric rectangular coordinates.
    x = (N+h) * cosf * cosl;
    y = (N+h) * cosf * sinl;
    z = ((1.0e0-e2) * N + h) * sinf;

    // Finished.
    return;
}

void
ell2car(double phi, double lambda, double h, double& x, double& y, double& z, 
    const Ellipsoid& e)
noexcept
{
    // Eccentricity squared.
    double e2 { e.eccentricity_squared() };

    // Radius of curvature in the prime vertical.
    double N { e.N(phi) };

    // Trigonometric numbers.
    double sinf { std::sin(phi) };
    double cosf { std::cos(phi) };
    double sinl { std::sin(lambda) };
    double cosl { std::cos(lambda) };

    // Compute geocentric rectangular coordinates.
    x = (N+h) * cosf * cosl;
    y = (N+h) * cosf * sinl;
    z = ((1.0e0-e2) * N + h) * sinf;

    // Finished.
    return;
}

} // end namespace

#endif
