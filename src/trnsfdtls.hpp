///
/// \file trnsfdtls.hpp
///
/// \brief Computational details for coordinate transformations.
///

#ifndef __NGPT_GEO_DETAILS__
#define __NGPT_GEO_DETAILS__

namespace ngpt
{

namespace detail
{

/// Compute the rotation matrix R, to convert a cartesian vector to a
/// topocentric one, i.e.
/// [n,e,u]**T = R*[dx, dy, dz]**T
///
/// \param[in] sinf sin(latitude)
/// \param[in] sinl sin(longtitude)
/// \param[in] cosf cos(latitude)
/// \param[in] cosl cos(longtitude)
/// \param[in] coef At output holds the coefficients of the rotation matrix R,
///                 in row major form; must have size >= 9.
/// \throw          Does not throw.
///
void
car2top_matrix(double sinf, double sinl, double cosf, double cosl, double* coef)
noexcept;

/// Compute the rotation matrix R, to convert a cartesian vector of variances
/// to topocentric one, i.e.
/// [sn**2,se**2,su**2]**T = R*[sx**2, sy**2, sz**2]**T
///
/// \param[in] sinf sin**2(latitude)
/// \param[in] sinl sin**2(longtitude)
/// \param[in] cosf cos**2(latitude)
/// \param[in] cosl cos**2(longtitude)
/// \param[in] coef At output holds the coefficients of the rotation matrix R,
///                 in row major form; must have size >= 9.
/// \throw          Does not throw.
///
void
car2top_cov_matrix(double sin2f, double sin2l, double cos2f, double cos2l,
    double* coef) noexcept;
} // namespace detail

} // namespace ngpt

#endif
