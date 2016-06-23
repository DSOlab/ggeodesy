///
/// \file car2ell.hpp
///
/// \brief Transformation of ellipsoidal to cartesian coordinates.
///

#ifndef __CARTESIAN_TO_ELLIPSOIDAL__
#define __CARTESIAN_TO_ELLIPSOIDAL__

#include <cmath>
#include "ellipsoid.hpp"
#include "geoconst.hpp"  // Needed for ngpt::DPI

namespace ngpt
{

/// \brief Cartesian to ellipsoidal.
///
/// Transform cartesian, geocentric coordinates (x, y, z) to ellipsoidal (i.e.
/// latitude, longtitude, ellispoidal height). This is a template function, 
/// depending on the ellipsoid parameters. All units are meters and radians.
///
/// \tparam     E      The reference ellipsoid (i.e. one of ngpt::ellipsoid).
/// \param[in]  x      Cartesian x-component, meters.
/// \param[in]  y      Cartesian y-component, meters.
/// \param[in]  z      Cartesian z-component, meters.
/// \param[out] phi    Ellipsoidal latitude, radians.
/// \param[out] lambda Ellipsoidal longtitude, radians.
/// \param[out] h      Ellipsoidal height, meters.
/// \throw             Does not throw.
///
/// \see Fukushima, T., "Transformation from Cartesian to geodetic coordinates
///      accelerated by Halley's method", J. Geodesy (2006), 79(12): 689-693
///
template<ellipsoid E>
    void car2ell(double x, double y, double z, double& phi, double& lambda,
        double& h) noexcept
{
    constexpr double a { ellipsoid_traits<E>::a };
    constexpr double f { ellipsoid_traits<E>::f };
    
    // Functions of ellipsoid parameters.
    constexpr double aeps2 { a*a*1e-32 };
    constexpr double e2    { (2.0e0-f)*f };
    constexpr double e4t   { e2*e2*1.5e0 };
    constexpr double ep2   { 1.0e0-e2 };
    constexpr double ep    { std::sqrt(ep2) };
    constexpr double aep   { a*ep };

    // Compute Coefficients of (Modified) Quartic Equation
    // Remark: Coefficients are rescaled by dividing by 'a'

    // Compute distance from polar axis squared.
    double p2 { x*x + y*y };

    // Compute longitude lambda.
    if ( p2 ) {
        lambda = std::atan2(y,x);
    } else { 
        lambda = .0e0;
    }

    // Ensure that Z-coordinate is unsigned.
    double absz { std::abs(z) };

    if (p2 > aeps2) { // Continue unless at the poles

        // Compute distance from polar axis.
        double p { std::sqrt(p2) };
        // Normalize.
        double s0  { absz/a };
        double pn  { p/a };
        double zp  { ep*s0 };
        // Prepare Newton correction factors.
        double c0  { ep*pn };
        double c02 { c0*c0 };
        double c03 { c02*c0 };
        double s02 { s0*s0 };
        double s03 { s02*s0 };
        double a02 { c02+s02 };
        double a0  { std::sqrt(a02) };
        double a03 { a02*a0 };
        double d0  { zp*a03 + e2*s03 };
        double f0  { pn*a03 - e2*c03 };
        // Prepare Halley correction factor.
        double b0  { e4t*s02*c02*pn*(a0-ep) };
        double s1  { d0*f0 - b0*s0 };
        double cp  { ep*(f0*f0-b0*c0) };
        // Evaluate latitude and height.
        phi = ::atan(s1/cp);
        double s12 { s1*s1 };
        double cp2 { cp*cp };
        h = (p*cp+absz*s1-a*std::sqrt(ep2*s12+cp2))
            /std::sqrt(s12+cp2);

      } else { // Special case: pole.

        phi = ngpt::DPI / 2e0;
        h   = absz - aep;
    }

    // Restore sign of latitude.
    if (z < 0.e0) phi = -phi;

    // Finished.
    return;
}

} // end namespace

#endif
