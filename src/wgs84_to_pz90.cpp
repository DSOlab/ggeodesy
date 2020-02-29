///
/// @file wgs84_to_pz90.cpp
///

#include "geodesy.hpp"

void
ngpt::wgs84_to_pz90(const double *xwgs, double *xpz, int pts, int selection)
{
  auto prms = ngpt::wgs2pz_parameters[selection];
  for (int i=0; i<pts; i++) {
    int ofst = i*3;
    xpz[ofst+0] = xwgs[ofst+0] + prms.tx + prms.d*xwgs[ofst+0] 
                - prms.r3*xwgs[ofst+1] + prms.r2*xwgs[ofst+2];
    xpz[ofst+1] = xwgs[ofst+1] + prms.ty + prms.r3*xwgs[ofst+0]
                + prms.d*xwgs[ofst+1]  - prms.r1*xwgs[ofst+2];
    xpz[ofst+2] = xwgs[ofst+2] + prms.tz - prms.r2*xwgs[ofst+0] 
                + prms.r1*xwgs[ofst+1] - prms.d*xwgs[ofst+2];
  }
  return;
}
