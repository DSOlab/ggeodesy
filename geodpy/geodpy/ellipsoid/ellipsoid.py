#! /usr/bin/python
from __future__ import print_function
import math

class Ellipsoid:
    ''' A class to represent a Reference Ellipsoid. It can be a "standard" ellipsoid,
      e.g. WGS84, GRS80, etc., or a user-defined one. This depends on the instance
      construction (see the constructor or whatever Pythons calls constructors by)
    '''

    __known_ell__ = {'grs80':{'a':6378137.0e0, 'f':1e0/298.257222101e0},
        'wgs84':{'a':6378137.0e0, 'f':1e0/298.257223563e0},
        'pz90' :{'a':6378135.0e0, 'f':1e0/298.257839303e0}}

    def __init__(self, *args, **kwargs):
        ''' Constructing an ellipsoid can be done in one of two ways :

        * Construct a standard ellipsoid, e.g. ::

            refel = Ellipsoid('grs80')

          use ``Ellipsoid('name')``, where ``'name'`` is one of:

          #. 'grs80'
          #. 'wgs84'
          #. 'pz90'

        * Construct a user-defined ellipsoid, e.g ::

            refel = Ellipsoid(a=6378137.0, f=1e0/298.257222101e0, name='my ell'),

          where the parameter ``name`` is optional.

        ''' 
        ## Construct from ellipsoid name
        if len(args) == 1:
            name = args[0].lower()
            if name not in self.__known_ell__:
                raise RuntimeError('Invalid ellipsoid name {:s}'.format(name))
            self.__a = self.__known_ell__[name]['a']
            self.__f = self.__known_ell__[name]['f']
            self.__n = name
        ## Construct user-defined ellipsoid
        elif len(kwargs) >=2 and len(kwargs) < 4:
            self.__a = kwargs.get('a', None)
            self.__f = kwargs.get('f', None)
            if self.__a is None or self.__f is None:
                raise RuntimeError('Invalid ellipsoid initialization')
            self.__n = kwargs.get('name', 'user-defined')
        else:
            raise RuntimeError('Invalid ellipsoid initialization. Read the docs')

    def eccentricity_squared(self):
        return (2e0 - self.__f) * self.__f

    def semi_major(self):
        return self.__a

    a = semi_major

    def flattening(self):
        return self.__f

    f = flattening

    def semi_minor(self):
        return self.__a * (1e0 - self.__f)

    b = semi_minor

    def N(self, lat_in_rad):
        ''' Compute the normal radius of curvature at a given latitude.
        '''
        cosf  = math.cos(lat_in_rad)
        sinf  = math.sin(lat_in_rad)
        acosf = self.__a * cosf
        bsinf = self.semi_minor() * sinf
        den   = math.sqrt(acosf*acosf + bsinf*bsinf)
        return (self.__a * self.__a) / den

    def M(self, lat_in_rad):
        ''' Compute the meridional radius of curvature at a given latitude.
        '''
        cosf  = math.cos(lat_in_rad)
        sinf  = math.sin(lat_in_rad)
        acosf = self.__a * cosf
        bsinf = self.semi_minor() * sinf
        tmpd  = acosf*acosf + bsinf*bsinf
        return ((a*b)/tmp) * ((self.__a*self.semi_minor())/math.sqrt(tmpd))

if __name__ == "__main__":
    ## Ellipsoid initialization
    el1 = Ellipsoid('wgs84')
    try:
        el2 = Ellipsoid('wgs48')
    except:
        print('Exception caught')
        pass
    el3 = Ellipsoid(a=1, f=2, name='bar')
