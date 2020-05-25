import math
import ellipsoid

def __haversine__(angle):
  sinHalfTheta = math.sin(angle/2.0e0)
  return sinHalfTheta * sinHalfTheta

def haversine(lat1, lon1, lat2, lon2, ell):
  h = __haversine__(lat2-lat1) + math.cos(lat1)*math.cos(lat2)*__haversine__(lon2-lon1)
  EarthRadius = (2e0*ell.a() + ell.semi_minor())/3e0
  return 2e0 * EarthRadius * math.asin(math.sqrt(h))


def inverse_vincenty(lat1, lon1, lat2, lon2, ell, convergence_limit=1e-10):
  MAX_ITERATIONS = 100;
  iteration = 0;

  a = ell.a()
  f = ell.f()
  b = ell.semi_minor()

  U1     = math.atan((1-f)*math.tan(lat1))
  U2     = math.atan((1-f)*math.tan(lat2))
  L      = lon2 -lon1
  sinU1  = math.sin(U1)
  sinU2  = math.sin(U2)
  cosU1  = math.cos(U1)
  cosU2  = math.cos(U2)
  lmbda  = L
  lmbdaP = float("inf")

  while (abs(lmbda-lmbdaP)) > convergence_limit:
    if iteration > MAX_ITERATIONS:
      raise RuntimeError("Inverse Vincenty cannot converge")
    sinLambda  = math.sin(lmbda)
    cosLambda  = math.cos(lmbda)
    sinSigma   = math.sqrt( (cosU2*sinLambda) * (cosU2*sinLambda) + 
      (cosU1*sinU2-sinU1*cosU2*cosLambda) * (cosU1*sinU2-sinU1*cosU2*cosLambda) )
    cosSigma   = sinU1*sinU2 + cosU1*cosU2*cosLambda
    sigma      = math.atan2(sinSigma, cosSigma)
    sinAlpha   = cosU1*cosU2*sinLambda / sinSigma
    cosSqAlpha = 1e0 - sinAlpha*sinAlpha
    cos2SigmaM = cosSigma - 2e0*sinU1*sinU2/cosSqAlpha
    C          = (f/16e0)*cosSqAlpha*(4e0+f*(4e0-3e0*cosSqAlpha))
    lmbdaP     = lmbda
    lmbda      = L + (1e0-C)*f*sinAlpha*(sigma+C*sinSigma*
      (cos2SigmaM+C*cosSigma*(-1e0+2e0*cos2SigmaM*cos2SigmaM)))

  uSq = cosSqAlpha*(a*a-b*b)/(b*b)
  k1  = (math.sqrt(1e0+uSq)-1e0)/(math.sqrt(1e0+uSq)+1e0)
  A   = (1e0+0.25e0*k1*k1)/(1e0-k1)
  B = k1*(1e0-(3e0/8e0)*k1*k1)
  deltaSigma = B*sinSigma*(cos2SigmaM+B/4e0*(cosSigma*
    (-1e0+2e0*cos2SigmaM*cos2SigmaM)-B/6e0*cos2SigmaM*(-3e0+4e0*sinSigma*sinSigma)
    *(-3e0+4e0*cos2SigmaM*cos2SigmaM)))
  distance   = b*A*(sigma-deltaSigma)
    
  # forward azimouth
  a12 = math.atan2(cosU2*sinLambda, cosU1*sinU2-sinU1*cosU2*cosLambda)
  # normalize
  a12 = math.fmod(a12+2e0*math.pi, 2*math.pi)

  # backward azimouth
  a21 = math.atan2(cosU1*sinLambda, -sinU1*cosU2+cosU1*sinU2*cosLambda)
  # normalize
  a21 = math.fmod(a21+2e0*math.pi, 2*math.pi)
    
  return distance, a12, a21
