#! /usr/bin/python
from __future__ import print_function
import math

from geodpy.ellipsoid.ellipsoid import Ellipsoid
import geodpy.crdtrans as gcrd
import geodpy.units.units as gunt

def dist(lst3): return math.sqrt(lst3[0]*lst3[0]+lst3[1]*lst3[1]+lst3[2]*lst3[2])
def sublst(lst1, lst2): return [ lst1[i]-lst2[i] for i in range(len(lst1)) ]

dyng = [4595220.009e0, 2039434.082e0, 3912626.005e0]
rand = [ i+5e0 for i in dyng ]

ell = Ellipsoid('grs80')

lat, lon, hgt = gcrd.cartesian2ellipsoidal(dyng[0], dyng[1], dyng[2], ell)
x, y, z       = gcrd.ellipsoidal2cartesian(lon, lat, hgt, ell)
dx, dy, dz    = gcrd.cartesian2topocentric(dyng[0], dyng[1], dyng[2], x=rand[0], y=rand[1], z=rand[2], ellipsoid=ell)

print('Coordinate transformations:')
print('Dyng ellipsoidal coordinates:')
d, m, s = gunt.rad2hexdeg(lat)
print('\tLatitude  : {:+3d} {:2d} {:8.5f}'.format(d, m, s))
d, m, s = gunt.rad2hexdeg(lon)
print('\tLongtitude: {:+3d} {:2d} {:8.5f}'.format(d, m, s))
print('\tHeight    : {:+10.4f}'.format(hgt))

print('Differences:')
print('\tDx={:10.5f}\n\tDy={:10.5f}\n\tDz={:10.5f}'.format(abs(dyng[0]-x), abs(dyng[1]-y), abs(dyng[2]-z)))
assert(abs(dyng[0]-x)<1e-12 and abs(dyng[1]-y)<1e-12 and abs(dyng[2]-z)<1e-12)

print('Topocentric Vector')
print('\tDn={:10.5f}\n\tDe={:10.5f}\n\tDu={:10.5f}\n\tS={:5.1f}'.format(dx, dy, dz, math.sqrt(dx*dx+dy*dy+dz*dz)))
assert(abs(dist([dx, dy, dz])-dist(sublst(rand,dyng)))<1e-12)
