///
/// \file vincenty.hpp
///
/// \brief Implementation of Vincenty's formulae and other algorithms to
///        compute great circle/godesic quantities.
///
/// \todo
///       - Need to test!
///       - Implement the direct/inverse algorithms from GographicLib, see
///         http://geographiclib.sourceforge.net/html/classGeographicLib_1_1Geodesic.html

#ifndef __NGPT_VINCENTY__HPP
#define __NGPT_VINCENTY_HPP

#include <cmath>
#include <stdexcept>
#include "ellipsoid.hpp"
#ifdef DEBUG
#include <iostream>
#include <iomanip>
#endif

namespace ngpt {

/// \brief Compute the haversine function.
///
/// This function is used within the haversine formula to compute great-circle
/// distances between two points on a sphere from their longitudes and
/// latitudes.
///
/// \see https://en.wikipedia.org/wiki/Haversine_formula
double haversine(double angle) noexcept
{
    double sinHalfTheta { std::sin(angle/2.0) };
    return sinHalfTheta * sinHalfTheta;
}

/// \brief Compute the great circle distance between two points using the
/// haversine formula.
///
/// Given the (ellipsoidal) coordinates of two points (1 and 2) in radians,
/// compute the great circle distance between them, using the haversine formula
/// see https://en.wikipedia.org/wiki/Haversine_formula.
///
/// For this computation we need the radius of the sphere/ellipsoid. Here, we
/// considere the mean radius as: 2a+b/3
///
/// \warning h only approaches 1 for antipodal points (on opposite sides of the
///          sphere). This will cause the formula to fail!
///
/// \note The haversine formula and law of cosines can't be guaranteed correct
///       to better than 0.5%
///
/// \see https://en.wikipedia.org/wiki/Haversine_formula
///
template<ellipsoid E = ellipsoid::wgs84>
    double haversine(double lat1, double lon1, double lat2, double lon2)
{
    double h  { haversine(lat2-lat1) + std::cos(lat1) * std::cos(lat2) * 
        haversine(lon2-lon1) };
    // according to wikipedia, The International Union of Geodesy and Geophysics
    // (IUGG) defines the mean radius (denoted R_1) to be
    // 2a+b/3
    double EarthRadius { (2*ellipsoid_traits<E>::a + semi_minor<E>())/3 };
    return 2 * EarthRadius * std::asin(std::sqrt(h));
}


/// \brief Compute the inverse Vincenty formula.
///
/// Given the (ellipsoidal) coordinates of two points (1 and 2) in radians,
/// calculate the forward azimouths from 1->2 (i.e. a12) and from 2->1 (i.e. a21)
/// and the ellipsoidal (along the geodesic) distance between the two points.
/// The inverse Vincenty's formula is used for the computation.
///
/// \throw Will throw an std::out_of_range (implying the input aruments/points)
///        in case more than MAX_ITERATIONS(=100) are needed for the algorithm
///        to converge.
///
/// \see https://en.wikipedia.org/wiki/Vincenty's_formulae
///
template<ellipsoid E = ellipsoid::wgs84>
    double inverse_vincenty(double lat1, double lon1, double lat2, double lon2,
        double& a12, double& a21, double convergence_limit = 1e-10)
{
    const int MAX_ITERATIONS = 100;
    int iteration = 0;

    double a = ellipsoid_traits<E>::a;
    double f = ellipsoid_traits<E>::f;
    double b = semi_minor<E>();

    double U1     { std::atan((1-f)*std::tan(lat1)) };
    double U2     { std::atan((1-f)*std::tan(lat2)) };
    double L      { lon2 -lon1 };
    double sinU1  { std::sin(U1) };
    double sinU2  { std::sin(U2) };
    double cosU1  { std::cos(U1) };
    double cosU2  { std::cos(U2) };
    double lambda {L};
    double sinSigma, cosSigma, sigma, sinAlpha, cosSqAlpha, C, cos2SigmaM,
           lambdaP, sinLambda, cosLambda;

    do {
        if (++iteration > MAX_ITERATIONS) {
            throw std::out_of_range("Inverse Vincenty cannot converge after 100 iterations!");
        }
        sinLambda  = std::sin(lambda);
        cosLambda  = std::cos(lambda);
        sinSigma   = std::sqrt( (cosU2*sinLambda) * (cosU2*sinLambda) +
                (cosU1*sinU2-sinU1*cosU2*cosLambda) * (cosU1*sinU2-sinU1*cosU2*cosLambda) );
        cosSigma   = sinU1*sinU2 + cosU1*cosU2*cosLambda;
        sigma      = atan2(sinSigma, cosSigma);
        sinAlpha   = cosU1*cosU2*sinLambda / sinSigma;
        cosSqAlpha = 1.0 - sinAlpha*sinAlpha;
        cos2SigmaM = cosSigma - 2.0*sinU1*sinU2/cosSqAlpha;
        C          = (f/16.0)*cosSqAlpha*(4.0+f*(4.0-3.0*cosSqAlpha));
        lambdaP    = lambda;
        lambda     = L + (1.0-C)*f*sinAlpha*
                    (sigma+C*sinSigma*(cos2SigmaM+C*cosSigma*
                    (-1.0+2.0*cos2SigmaM*cos2SigmaM)));
    } while ( std::abs(lambda-lambdaP) > convergence_limit );
#ifdef DEBUG
    std::cout<<"\t[DEBUG]Vincenty Inverse converged after "<<iteration<<" iterations\n";
#endif

    double uSq { cosSqAlpha*(a*a-b*b)/(b*b) };
    double k1  { (std::sqrt(1.0+uSq)-1.0)/(std::sqrt(1.0+uSq)+1.0) };
    double A   { (1+0.25*k1*k1)/(1.0-k1) };
    double B   { k1*(1.0-(3.0/8.0)*k1*k1) };
    double deltaSigma { B*sinSigma*(cos2SigmaM+B/4.0*(cosSigma*
            (-1.0+2.0*cos2SigmaM*cos2SigmaM)-B/6.0*cos2SigmaM*
            (-3.0+4.0*sinSigma*sinSigma)*(-3.0+4.0*cos2SigmaM*cos2SigmaM))) };
    double distance { b*A*(sigma-deltaSigma) };
    
    // forward azimouth
    a12 = std::atan2(cosU2*sinLambda, cosU1*sinU2-sinU1*cosU2*cosLambda);
    // normalize
    a12 = std::fmod(a12+D2PI, D2PI);

    // backward azimouth
    a21 = std::atan2(cosU1*sinLambda, -sinU1*cosU2+cosU1*sinU2*cosLambda);
    // normalize
    a21 = std::fmod(a21+DPI, D2PI);
    
    return distance;
}

/// \brief Direct Vincenty formula.
///
/// Given an initial point (lat1, lon1) and initial azimuth, a1, and a distance,
/// s, along the geodesic the problem is to find the end point (lat2, lon2)
/// and azimuth, a2.
///
/// \warning 
///
/// \see https://en.wikipedia.org/wiki/Vincenty%27s_formulae
template<ellipsoid E = ellipsoid::wgs84>
    double direct_vincenty(double lat1, double lon1, double a1, double s, 
        double& lat2, double& lon2, double convergence_limit = 1e-10)
{
    const int MAX_ITERATIONS = 100;
    int iteration = 0;

    double a = ellipsoid_traits<E>::a;
    double f = ellipsoid_traits<E>::f;
    double b = semi_minor<E>();
    
    double cosa1  { std::cos(a1) };
    double sina1  { std::sin(a1) };
    double U1     { std::atan((1-f)*std::tan(lat1)) };
    double sigma1 { std::atan2(std::tan(U1), cosa1) };
    double sina   { std::cos(U1) * sina1 };
    double cosaSq {1.0e0-sina*sina};
    double uSq    { cosaSq*(a*a-b*b)/(b*b) };
    double k1     { (std::sqrt(1.0+uSq)-1.0)/(std::sqrt(1.0+uSq)+1.0) };
    double A      { (1.0+0.25*k1*k1)/(1.0-k1) };
    double B      { k1*(1.0-(3.0/8.0)*k1*k1) };

    double sigma { s/b*A }; // initial guess
    double sigmaP, sigmaM2, cosSigmaM2, deltaSigma, sinSigma;

    do {
        if (++iteration > MAX_ITERATIONS) {
            throw std::out_of_range("Direct Vincenty cannot converge after 100 iterations!");
        }
        sigmaM2    = 2.0*sigma1 + sigma;
        cosSigmaM2 = std::cos(sigmaM2);
        sinSigma   = std::sin(sigma);
        deltaSigma = B*sinSigma*(cosSigmaM2 
            + (1/4.0)*B*(std::cos(sigma)*(-1+2*cosSigmaM2*cosSigmaM2)
            - (1/6.0)*B*cosSigmaM2*(-3+4*sinSigma*sinSigma)*(-3+4*cosSigmaM2*cosSigmaM2)));
        sigmaP = sigma;
        sigma = s/(b*A) + deltaSigma;
    } while ( std::abs(sigma-sigmaP) > convergence_limit );
#ifdef DEBUG
    std::cout<<"\t[DEBUG]Direct Vincenty converged after "<<iteration<<" iterations\n";
#endif

    double sinU1    { std::sin(U1) };
    double cosU1    { std::cos(U1) };
    double cosSigma { std::cos(sigma) };

    // compute latitude
    double nom   { sinU1*cosSigma+cosU1*sinSigma*cosa1 };
    double denom { sina*sina+std::pow(sinU1*sinSigma-cosU1*cosSigma*cosa1, 2.0) };
    denom        = std::sqrt(denom)*(1.0-f);
    lat2         = std::atan2(nom, denom);

    // compute longtitude
    nom           = sinSigma*sina1;
    denom         = cosU1*cosSigma - sinU1*sinSigma*cosa1;
    double lambda { std::atan2(nom, denom) };
    double C      { (f/16.0)*cosaSq*(4.0+f*(4.0-3.0*cosaSq)) };
    double L      { lambda - (1.0-C)*f*sina*(sigma+C*sinSigma*(cosSigmaM2+C*cosSigma*
                    (-1.0+2.0*cosSigmaM2*cosSigmaM2))) };
    lon2          = L + lon1;
    
    // compute azimouth
    return std::atan2(sina, -sinU1*sinSigma+cosU1*cosSigma*cosa1) + DPI;
}

} // end namespace ngpt
#endif
