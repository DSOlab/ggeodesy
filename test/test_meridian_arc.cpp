#include "ellipsoid.hpp"
#include "geodesy.hpp"
#include <iostream>

int main() {
  double lat = ngpt::deg2rad(60e0);
  
  double s = ngpt::meridian_arc_length<ngpt::ellipsoid::grs80>(lat,0);
  printf("\nMeridian Arc Length to lat = %+5.1f is %+20.10fm", ngpt::rad2deg(lat), s);
  
  s = ngpt::meridian_arc_length<ngpt::ellipsoid::grs80>(lat,1);
  printf("\nMeridian Arc Length to lat = %+5.1f is %+20.10fm", ngpt::rad2deg(lat), s);

  s = ngpt::meridian_arc_length<ngpt::ellipsoid::grs80>(lat,2);
  printf("\nMeridian Arc Length to lat = %+5.1f is %+20.10fm", ngpt::rad2deg(lat), s);
  
  printf("\n");
  return 0;
}
