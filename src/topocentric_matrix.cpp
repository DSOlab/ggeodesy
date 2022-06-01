#include "geodesy.hpp"

using dso::MATRIX3x3;

/// @brief Given the geodetic coordinates of a reference point, return the
///        matrix that turns any vector from the ference point to point P to 
///        topocentric coordinates (e,n,u)
MATRIX3x3 dso::topocentric_matrix(double lambda, double phi) noexcept {
  const double cf = std::cos(phi);
  const double sf = std::sin(phi);
  const double cl = std::cos(lambda);
  const double sl = std::sin(lambda);
  
  // unit vector along east
  const double d00 = -sl;
  const double d01 = cl;
  const double d02 = 0e0;
  // unit vector along north
  const double d10 = -sf * cl;
  const double d11 = -sf * sl;
  const double d12 = cf;
  // unit vector along up
  const double d20 = cf*cl;
  const double d21 = cf*sl;
  const double d22 = sf;

#ifdef USE_EIGEN
  double data[9] = {d00,d01,d02,d10,d11,d12,d20,d21,d22};
  return Eigen::Map<MATRIX3x3>(data,3,3);
#else
  return MATRIX3x3({d00,d01,d02,d10,d11,d12,d20,d21,d22});
#endif
}
