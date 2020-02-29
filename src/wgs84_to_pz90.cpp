///
/// @file wgs84_to_pz90.cpp
///

#include "geodesy.hpp"

void
ngpt::pz90_to_wgs84(const double *xwgs, double *xpz, int pts, int selection)
{
  auto prms = ngpt::pz2wgs_parameters[selection];
  // angles to rad*1e9
  double r1 = prms.r1*4.847309743e0;
  double r2 = prms.r2*4.847309743e0;
  double r3 = prms.r3*4.847309743e0;
  // double coordinates in meters*1e-9
  double xMm, yMm, zMm;
  for (int i=0; i<pts; i++) {
    int ofst = i*3;
    xMm = xwgs[ofst+0] * 1e-9;
    yMm = xwgs[ofst+1] * 1e-9;
    zMm = xwgs[ofst+2] * 1e-9;
    xpz[ofst+0] = xwgs[ofst+0] + prms.tx + prms.d*xMm - r3*yMm + r2*zMm;
    xpz[ofst+1] = xwgs[ofst+1] + prms.ty + r3*xMm + prms.d*yMm  - r1*zMm;
    xpz[ofst+2] = xwgs[ofst+2] + prms.tz - r2*xMm + r1*yMm - prms.d*zMm;
  }
  return;
}
