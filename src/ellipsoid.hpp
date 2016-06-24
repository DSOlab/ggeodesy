///
/// \file ellipsoid.hpp
///
/// \brief Definition and basic setup of an ellipsoid (kinda) class.
///
/// This file defines a list of frequently used (referennce) Ellipsoids in
/// geodesy. Each ellipsoid comes with a list of fundamental geometric
/// characteristics (semi-major axis, flattening, name).
///
/// The ellipsoids are not a class; rather an enumeration (ngpt::ellipsoid).
/// A (per-ellipsoid) specialized traits class (i.e. ngpt::ellipsoid_traits)
/// is used to implement the basic properties of each ellipsoid.
///
/// To add a new ellipsoid, add a new item in the ngpt::ellipsoid enum, and
/// provide its characteristics via a specialization of the ngpt::ellipsoid_traits
/// class
///
/// \example test_geodesy.cpp
///

#ifndef __REFERENCE_ELLIPSOID__
#define __REFERENCE_ELLIPSOID__

#include <cmath>
#ifdef DEBUG
#include <iostream>
#endif

namespace ngpt
{

/// \brief A list of well-known reference ellipsoids.
///
/// For each of these reference ellipsoids, a series of traits (i.e. geometric
/// characteristics) will be specialized later on, using the template class
/// ngpt::ellipsoid_traits.
///
enum class ellipsoid : char
{
    grs80,
    wgs84,
    pz90
};

/// \brief A class to hold ellipsoid traits (generic case).
///
/// A (template class) to hold specialized geometric quantities for each
/// of the eumerated elements (i.e. reference ellipsoids) in the 
/// ngpt::ellipsoid enum.
/// I.e., to make any element of ngpt::ellipsoid usable, specialize this
/// (trait) class.
///
/// \tparam E  The reference ellipsoid to be specialized (i.e. one of 
///            ngpt::ellipsoid).
template<ellipsoid E> struct ellipsoid_traits {};

/// \brief A class to hold traits for the GRS-80 (i.e. ngpt::ellispoid::grs80)
/// reference ellipsoid.
///
/// \see https://en.wikipedia.org/wiki/GRS_80
template<>
    struct ellipsoid_traits<ellipsoid::grs80>
{
    /// Semi-major axis (m).
    static constexpr double a      { 6378137.0e0 };
    /// Flattening.
    static constexpr double f      { 1.0e00/298.257222101e0 };
    /// Reference ellipsoid name.
    static constexpr const char* n { "GRS80" };
};

/// \brief A class to hold traits for the WGS-84 (i.e. ngpt::ellispoid::wgs84)
/// reference ellipsoid.
///
/// \see https://en.wikipedia.org/wiki/World_Geodetic_System
template<>
    struct ellipsoid_traits<ellipsoid::wgs84>
{
    static constexpr double a      { 6378137.0e0 };
    static constexpr double f      { 1.0e00/298.257223563e0 };
    static constexpr const char* n { "WGS84" };
};

/// \brief A class to hold traits for the PZ-90 (i.e. ngpt::ellispoid::pz90)
/// reference ellipsoid.
///
/// \see http://www.navipedia.net/index.php/Reference_Frames_in_GNSS#GLONASS_reference_frame_PZ-90
template<>
    struct ellipsoid_traits<ellipsoid::pz90>
{
    static constexpr double a      { 6378135.0e0 };
    static constexpr double f      { 1.0e00/298.257839303e0 };
    static constexpr const char* n { "PZ90" };
};

/// \brief Compute the squared eccentricity.
///
/// For any reference ellipsoid (i.e. ngpt::ellipsoid) compute the squared 
/// (first) eccentricity (i.e. e^2).
///
/// \tparam E  The reference ellipsoid (i.e. one of ngpt::ellipsoid).
/// \return    Eccentricity squared.
/// \throw     Does not throw.
///
/// \see https://en.wikipedia.org/wiki/Geodetic_datum
template<ellipsoid E>
    constexpr double eccentricity_squared() noexcept
{
    constexpr double f { ellipsoid_traits<E>::f };
    return ( 2.0e0 - f ) * f;
}

/// \brief Compute the semi-minor axis (b).
///
/// For any reference ellipsoid (i.e. ngpt::ellipsoid) compute the semi-minor
/// axis 'b' in meters.
///
/// \tparam E The reference ellipsoid (i.e. one of ngpt::ellipsoid).
/// \return   Semi-minor axis of the reference ellipsoid in meters.
/// \throw    Does not throw.
///
/// \see https://en.wikipedia.org/wiki/Geodetic_datum
template<ellipsoid E>
    constexpr double semi_minor() noexcept
{
  constexpr double a { ellipsoid_traits<E>::a };
  constexpr double f { ellipsoid_traits<E>::f };
  return a * ( 1.0e0 - f );
}

/// \brief Compute the normal radius of curvature at a given latitude (on a 
///        reference ellipsoid).
///
/// \param[in] lat The latitude in radians.
/// \tparam    E   The reference ellipsoid (i.e. one of ngpt::ellipsoid).
/// \return        The normal radius of curvature in meters.
/// \throw         Does not throw (see notes).
///
/// \note If the denominator (den) is zero then funny things could happen; this
///       however should **never** occur for any reference ellipsoid.
///
/// \see "Physical Geodesy", pg. 194
/// \see https://en.wikipedia.org/wiki/Earth_radius
///
template<ellipsoid E>
    double N(double lat) noexcept
{
    constexpr double a { ellipsoid_traits<E>::a };
    double cosf  { std::cos(lat) };
    double sinf  { std::sin(lat) };
    double acosf { a * cosf };
    double bsinf { semi_minor<E>() * sinf };
    double den   { std::sqrt(acosf*acosf + bsinf*bsinf) };
#ifdef DEBUG
    double answer = (a*a)/den;
    double alternative = a/(std::sqrt(1.0-eccentricity_squared<E>()*sinf*sinf));
    if ( std::abs(answer-alternative) > 1e-8 ) {
        std::cerr<<"\n[ERROR] Too big a difference between normal radius of "
        << "curvature computation\nSee template<ellipsoid E> double N(double lat)"
        << "in file: ellipsoid.hpp\n";
        assert( false );
    }
#endif

    return (a * a) / den;
}

/// \brief Compute the meridional radii of curvature at a given latitude (on a 
///        reference ellipsoid).
///
/// \param[in] lat The latitude in radians.
/// \tparam    E   The reference ellipsoid (i.e. one of ngpt::ellipsoid).
/// \return        The meridional radius of curvature in meters.
/// \throw         Does not throw (see notes).
///
/// \see https://en.wikipedia.org/wiki/Earth_radius
///
template<ellipsoid E>
    double M(double lat) noexcept
{
    constexpr double a { ellipsoid_traits<E>::a };
    constexpr double b { semi_minor<E>() };
    double cosf  { std::cos(lat) };
    double sinf  { std::sin(lat) };
    double acosf { a * cosf };
    double bsinf { b * sinf };
    double tmpd  { acosf*acosf + bsinf*bsinf };
    return ( (a*b)/tmpd ) * ( (a*b)/std::sqrt(tmpd) );
}

} // end namespace ngpt

#endif
