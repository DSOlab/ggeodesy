#include <cmath>
#include "trnsfdtls.hpp"

/// @brief Rotation matrix for cartesian to topocentric transformation.
void
ngpt::detail::car2top_matrix(double sinf, double sinl, double cosf, double cosl,
    double* coef)
noexcept
{
    coef[0] = -sinf*cosl; coef[1] = -sinf*sinl; coef[2] = cosf;
    coef[3] = -sinl;      coef[4] = cosl;       coef[5] = 0e0;
    coef[6] = cosf*cosl;  coef[7] = cosf*sinl;  coef[8] = sinf;
}

/// @brief Rotation matrix for cartesian to topocentric transformation of
///        variances.
void
ngpt::detail::car2top_cov_matrix(double sin2f, double sin2l, double cos2f,
    double cos2l, double* coef)
noexcept
{
    coef[0] = sin2f*cos2l; coef[1] = sin2f*sin2l;  coef[2] = cos2f;
    coef[3] = sin2l;       coef[4] = cos2l;        coef[5] = 0e0;
    coef[6] = cos2f*cos2l; coef[7] = cos2f*sin2l;  coef[8] = sin2f;
}
