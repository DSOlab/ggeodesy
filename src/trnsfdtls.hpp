/// @file trnsfdtls.hpp
/// @brief Computational details for coordinate transformations.

#ifndef __NGPT_TRNSFDTLS_HPP__
#define __NGPT_TRNSFDTLS_HPP__

namespace dso {

namespace core {

/// Compute the rotation matrix R, to convert a cartesian vector to a
/// topocentric one, i.e.
/// \f$ [north, east, up]^T = R * [\delta x, \delta y, \delta z]^T \f$
/// The input trigonometric numbers (aka sin and cos of longtitude and
/// latitude) are computed using the ellipsoidal coordinates of the reference
/// point of the (resulting) topocentric system.
///
/// @param[in] sinf sin(latitude)
/// @param[in] sinl sin(longtitude)
/// @param[in] cosf cos(latitude)
/// @param[in] cosl cos(longtitude)
/// @param[in] coef At output holds the coefficients of the rotation matrix R,
///                 in row major form; must have size >= 9.
/// @throw          Does not throw.
///
void car2top_matrix(double sinf, double sinl, double cosf, double cosl,
                    double *coef) noexcept;

/// Compute the rotation matrix R, to convert a cartesian vector of variances
/// to topocentric one, i.e.
/// \f$ [Var_{north}, Var_{east}, Var_{up}]^T = R * [Var_{\delta x}, Var_{\delta
/// y}, Var_{\delta z}]^T \f$ The input trigonometric numbers (aka sin**2 and
/// cos**2 of longtitude and latitude) are computed using the ellipsoidal
/// coordinates of the reference point of the (resulting) topocentric system.
///
/// @param[in] sin2f sin**2(latitude) aka \f$ sin^2 \phi \f$
/// @param[in] sin2l sin**2(longtitude) aka \f$ sin^2 \lambda \f$
/// @param[in] cos2f cos**2(latitude) aka \f$ sin^2 \phi \f$
/// @param[in] cos2l cos**2(longtitude) aka \f$ cos^2 \lambda \f$
/// @param[in] coef  At output holds the coefficients of the rotation matrix R,
///                  in row major form; must have size >= 9.
/// @throw           Does not throw.
///
void car2top_cov_matrix(double sin2f, double sin2l, double cos2f, double cos2l,
                        double *coef) noexcept;
} // namespace core

} // namespace dso

#endif
