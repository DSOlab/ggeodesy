///
/// \file car2top.hpp
///
/// \brief Transformation of cartesian to topocentric vector.
///

#ifndef _CARTESIAN_TO_TOPOCENTRIC_
#define _CARTESIAN_TO_TOPOCENTRIC_

#include <cmath>

#include "car2ell.hpp"

namespace ngpt
{

/// \brief Cartesian to topocentric (vector).
///
/// Transform a vector expressed in cartesian coordinates to the topocentric,
/// local system around point i (i.e. North(i), East(i), Up(i)). This is a 
/// template function, depending on the ellipsoid parameter. All units in 
/// meters.
///
/// \tparam     E      The reference ellipsoid (i.e. one of ngpt::ellipsoid).
/// \param[in]  xi     Cartesian x-component of point i, meters.
/// \param[in]  yi     Cartesian y-component of point i, meters.
/// \param[in]  zi     Cartesian z-component of point i, meters.
/// \param[in]  xj     Cartesian x-component of point j, meters.
/// \param[in]  yj     Cartesian y-component of point j, meters.
/// \param[in]  zj     Cartesian z-component of point j, meters.
/// \param[out] north  Vector north component, meters.
/// \param[out] east   Vector east component, meters.
/// \param[out] up     Vector up component, meters.
/// \throw             Does not throw.
///
/// \note The ellispoid is needed to transform the cartesian coordinates of
///       the (reference) point i to ellispoidal coordinates.
///
/// \see "Physical Geodesy", pg. 209
///
template<ellipsoid E>
    void car2top(double xi, double yi, double zi, double xj, double yj, double zj, 
        double& north, double& east, double& up) noexcept
{

    // Ellipsoidal coordinates of reference point.
    double phi_i, lambda_i, h_i;

    // Cartesian to ellipsoidal for reference point.
    ngpt::car2ell<E>(xi, yi, zi, phi_i, lambda_i, h_i);

    // Trigonometric numbers.
    double cosf { std::cos(phi_i) };
    double cosl { std::cos(lambda_i) };
    double sinf { std::sin(phi_i) };
    double sinl { std::sin(lambda_i) };

    // Catresian vector.
    double dx { xj - xi };
    double dy { yj - yi };
    double dz { zj - zi };

    // Topocentric vector.
    north = - sinf * cosl * dx - sinf * sinl * dy + cosf * dz;
    east  = - sinl * dx        + cosl * dy;
    up    =   cosf * cosl * dx + cosf * sinl * dy + sinf * dz;

    // Finished.
    return;
}

} // end namespace

#endif
