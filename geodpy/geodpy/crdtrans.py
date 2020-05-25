from .ellipsoid.ellipsoid import Ellipsoid
import math

def ellipsoidal2cartesian(lon, lat, hgt, ellpsd=None):
    ''' Given a set of geocentric, ellipsoidal coordinates and optionaly a reference
      ellispoid, transform the set to cartesian coordinates, i.e. x, y, z.
    '''

    if ellpsd == None: ellpsd = Ellipsoid('grs80')

    ## Eccentricity squared.
    e2 = ellpsd.eccentricity_squared()

    ## Trigonometric numbers.
    sinf = math.sin(lat)
    cosf = math.cos(lat)
    sinl = math.sin(lon)
    cosl = math.cos(lon)

    ## Radius of curvature in the prime vertical.
    N = ellpsd.N(lat)

    ## Compute geocentric rectangular coordinates.
    x = (N+hgt) * cosf * cosl
    y = (N+hgt) * cosf * sinl
    z = ((1.0e0-e2) * N + hgt) * sinf

    return x, y, z

def cartesian2ellipsoidal(x, y, z, ellpsd=None):
    ''' Given a set of geocentric, cartesian coordinates and optionaly a reference
      ellispoid, transform the set to ellipsoidal coordinates, i.e. longtitude,
      latitude and (ellipsoidal) height.
    '''

    if ellpsd == None: ellpsd = Ellipsoid('grs80')

    ## Functions of ellipsoid parameters.
    a     = ellpsd.a()
    f     = ellpsd.f()
    aeps2 = a*a*1e-32
    e2    = ellpsd.eccentricity_squared()
    e4t   = e2*e2*1.5e0
    ep2   = 1.0e0-e2
    ep    = math.sqrt(ep2)
    aep   = a*ep

    ''' Compute Coefficients of (Modified) Quartic Equation
      Remark: Coefficients are rescaled by dividing by 'a'
    '''

    ## Compute distance from polar axis squared.
    p2 = x*x + y*y

    ## Compute longitude lon
    lon = math.atan2(y,x) if p2!=0e0 else 0e0

    ## Ensure that Z-coordinate is unsigned.
    absz = abs(z)

    if p2 > aeps2: ## Continue unless at the poles
        ## Compute distance from polar axis.
        p = math.sqrt(p2)
        ## Normalize.
        s0 = absz/a
        pn = p/a
        zp = ep*s0 
        ## Prepare Newton correction factors.
        c0  = ep*pn 
        c02 = c0*c0 
        c03 = c02*c0 
        s02 = s0*s0 
        s03 = s02*s0 
        a02 = c02+s02 
        a0  = math.sqrt(a02) 
        a03 = a02*a0 
        d0  = zp*a03 + e2*s03 
        f0  = pn*a03 - e2*c03 
        ## Prepare Halley correction factor.
        b0 = e4t*s02*c02*pn*(a0-ep) 
        s1 = d0*f0 - b0*s0 
        cp = ep*(f0*f0-b0*c0) 
        ## Evaluate latitude and height.
        lat = math.atan(s1/cp)
        s12 = s1*s1 
        cp2 = cp*cp 
        hgt = (p*cp+absz*s1-a*math.sqrt(ep2*s12+cp2))/math.sqrt(s12+cp2);
    else: ## Special case: pole.
        lat = math.pi / 2e0;
        hgt = absz - aep;

    ## Restore sign of latitude.
    if z < 0.e0: lat = -lat

    return lat, lon, hgt

def d_car2top(x, y, z, dx, dy, dz, ell):
    lat, lon, hgt = cartesian2ellipsoidal(x, y, z, ell)
    cosf = math.cos(lat)
    sinf = math.sin(lat)
    cosl = math.cos(lon)
    sinl = math.sin(lon)
    north = -sinf * cosl * dx - sinf * sinl * dy + cosf * dz;
    east  = -sinl * dx        + cosl * dy;
    up    =  cosf * cosl * dx + cosf * sinl * dy + sinf * dz;

    return north, east, up

def i_car2top(xi, yi, zi, xj, yj, zj, ell):
    # Cartesian vector.
    dx = xj - xi
    dy = yj - yi
    dz = zj - zi

    # transform to topocentric
    return d_car2top(xi, yi, zi, dx, dy, dz, ell)

def cartesian2topocentric(xref, yref, zref, **kwargs):
    if 'ellipsoid' in kwargs:
        ell = kwargs['ellipsoid']
    else:
        ell = Ellipsoid('grs80')
    if all(prs in kwargs for prs in ['dx', 'dy', 'dz']):
        return d_car2top(xref, yref, zref, kwargs['dx'], kwargs['dy'], kwargs['dz'], ell)
    if all(prs in kwargs for prs in ['x', 'y', 'z']):
        return i_car2top(xref, yref, zref, kwargs['x'], kwargs['y'], kwargs['z'], ell)
    raise RuntimeError('Invalid call to geodesy.car2top')
