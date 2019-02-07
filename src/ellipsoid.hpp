///
/// @file ellipsoid.hpp
///
/// @brief Definition and basic setup of an ellipsoid (kinda) class.
///
/// This file defines a list of frequently used (referennce) Ellipsoids in
/// geodesy. Each ellipsoid comes with a list of fundamental geometric
/// characteristics (semi-major axis, flattening, name).
///
/// There are two ways a user can make use of the ellipsoids:
/// * 1. If the ellipsoid of choice is known at compile time, then the
///      template function/traits can be used. A (per-ellipsoid) specialized
///      traits class (i.e. ngpt::ellipsoid_traits) is used to implement the
///      basic properties of each ellipsoid.
///      E.g., it the ellipsoid of choice is GRS80, then:
///      double semi_major = ellipsoid_traits<ellipsoid::grs80>::a;
///      or
///      double N = N<ellipsoid::grs80>(lat);
///      gives the normal radius of curvature at latitude lat.
///
/// * 2. If the ellipsoid of choice is only kown at runtime, then all relevant
///      computations/constants can be accessed via the Ellipsoid class; i.e.
///      Ellipsoid e (ellipsoid::grs80);
///      double semi_major = e.semi_major();
///      double N = e.N(lat);
///
/// @example test_geodesy.cpp
///

#ifndef __REFERENCE_ELLIPSOID__
#define __REFERENCE_ELLIPSOID__

#include <cmath>

namespace ngpt
{

/// @brief core namespace holds the core of ellipsoid-related functions
namespace core
{
    /// @brief Compute the squared eccentricity.
    ///
    /// Compute the squared (first) eccentricity (i.e. e^2) given the
    /// flattening of an ellipsoid, aka \f$ e^2 = (2-f)*f \f$
    /// @param[in] f flattening
    /// @return      squared eccentricity
    constexpr double
    eccentricity_squared(double f)
    noexcept
    { return ( 2.0e0 - f ) * f; }
    
    /// @brief Compute the semi-minor axis (aka b)
    ///
    /// Compute the semi-minor axis of an ellipsoid(i.e. e^2), given the
    /// flattening and the semi-major axis, aka \f$ \beta = \alpha * (1-f) \f$
    /// @param[in] f flattening
    /// @param[in] a semi-major axis (meters)
    /// @return      semi-minor axis (meters)
    constexpr double
    semi_minor(double a, double f)
    noexcept
    { return a * ( 1.0e0 - f ); }

    /// @brief Compute the normal radius of curvature at a given latitude (on 
    ///        a reference ellipsoid).
    ///
    /// @param[in] lat The latitude (radians)
    /// @param[in] a   The ellipsoid's semi-major axis (meters)
    /// @param[in] b   The ellipsoid's semi-minor axis (meters)
    /// @return        The normal radius of curvature (meters)
    /// @throw         Does not throw (see notes)
    ///
    /// @note If the denominator (den) is zero then funny things could happen;
    ///       this however should **never** occur for any reference ellipsoid.
    ///
    /// @see "Physical Geodesy", pg. 194
    /// @see https://en.wikipedia.org/wiki/Earth_radius
    ///
    constexpr double
    N(double lat, double a, double b) noexcept
    {
        double cosf  { std::cos(lat) };
        double sinf  { std::sin(lat) };
        double acosf { a * cosf };
        double bsinf { b * sinf };
        double den   { std::sqrt(acosf*acosf + bsinf*bsinf) };
        /*
        double answer = (a*a)/den;
        double alternative = a/(std::sqrt(1.0-eccentricity_squared<E>()*sinf*sinf));
        if ( std::abs(answer-alternative) > 1e-8 ) {
            std::cerr<<"\n[ERROR] Too big a difference between normal radius of "
            << "curvature computation\nSee template<ellipsoid E> double N(double lat)"
            << "in file: ellipsoid.hpp\n";
        }*/
        return (a * a) / den;
    }

    /// @brief Compute the meridional radii of curvature at a given latitude 
    ///        on a reference ellipsoid).
    ///
    /// @param[in] lat The latitude in (radians)
    /// @param[in] a   The ellipsoid's semi-major axis (meters)
    /// @param[in] b   The ellipsoid's semi-minor axis (meters)
    /// @return        The meridional radius of curvature (meters)
    /// @throw         Does not throw (see notes)
    ///
    /// @see https://en.wikipedia.org/wiki/Earth_radius
    ///
    constexpr double
    M(double lat, double a, double b) noexcept
    {
        double cosf  { std::cos(lat) };
        double sinf  { std::sin(lat) };
        double acosf { a * cosf };
        double bsinf { b * sinf };
        double tmpd  { acosf*acosf + bsinf*bsinf };
        return ( (a*b)/tmpd ) * ( (a*b)/std::sqrt(tmpd) );
    }
}// namespace core

/// @brief A list of well-known reference ellipsoids.
///
/// For each of these reference ellipsoids, a series of traits (i.e. geometric
/// characteristics) will be specialized later on, using the template class
/// ngpt::ellipsoid_traits.
///
enum class ellipsoid : char { grs80, wgs84, pz90 };

/// @brief A class to hold ellipsoid traits (generic case).
///
/// A (template class) to hold specialized geometric quantities for each
/// of the eumerated elements (i.e. reference ellipsoids) in the 
/// ngpt::ellipsoid enum.
/// I.e., to make any element of ngpt::ellipsoid usable, specialize this
/// (trait) class.
///
/// @tparam E  The reference ellipsoid to be specialized (i.e. one of 
///            ngpt::ellipsoid).
template<ellipsoid E> struct ellipsoid_traits {};

/// @brief A class to hold traits for the GRS-80 (i.e. ngpt::ellispoid::grs80)
/// reference ellipsoid.
///
/// @see https://en.wikipedia.org/wiki/GRS_80
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

/// @brief A class to hold traits for the WGS-84 (i.e. ngpt::ellispoid::wgs84)
/// reference ellipsoid.
///
/// @see https://en.wikipedia.org/wiki/World_Geodetic_System
template<>
    struct ellipsoid_traits<ellipsoid::wgs84>
{
    static constexpr double a      { 6378137.0e0 };
    static constexpr double f      { 1.0e00/298.257223563e0 };
    static constexpr const char* n { "WGS84" };
};

/// @brief A class to hold traits for the PZ-90 (i.e. ngpt::ellispoid::pz90)
/// reference ellipsoid.
///
/// @see http://www.navipedia.net/index.php/Reference_Frames_in_GNSS#GLONASS_reference_frame_PZ-90
template<>
    struct ellipsoid_traits<ellipsoid::pz90>
{
    static constexpr double a      { 6378135.0e0 };
    static constexpr double f      { 1.0e00/298.257839303e0 };
    static constexpr const char* n { "PZ90" };
};

/// @brief Compute the squared eccentricity.
///
/// For any reference ellipsoid (i.e. ngpt::ellipsoid) compute the squared 
/// (first) eccentricity (i.e. e^2).
///
/// @tparam E  The reference ellipsoid (i.e. one of ngpt::ellipsoid).
/// @return    Eccentricity squared.
/// @throw     Does not throw.
///
/// @see https://en.wikipedia.org/wiki/Geodetic_datum
template<ellipsoid E>
    constexpr double
    eccentricity_squared() noexcept
{
    return core::eccentricity_squared(ellipsoid_traits<E>::f);
}

/// @brief Compute the semi-minor axis (b).
///
/// For any reference ellipsoid (i.e. ngpt::ellipsoid) compute the semi-minor
/// axis 'b' in meters.
///
/// @tparam E The reference ellipsoid (i.e. one of ngpt::ellipsoid).
/// @return   Semi-minor axis of the reference ellipsoid in meters.
/// @throw    Does not throw.
///
/// @see https://en.wikipedia.org/wiki/Geodetic_datum
template<ellipsoid E>
    constexpr double
    semi_minor() noexcept
{
  return core::semi_minor(ellipsoid_traits<E>::a, ellipsoid_traits<E>::f);
}

/// @brief Compute the normal radius of curvature at a given latitude (on a 
///        reference ellipsoid).
///
/// @param[in] lat The latitude in radians.
/// @tparam    E   The reference ellipsoid (i.e. one of ngpt::ellipsoid).
/// @return        The normal radius of curvature in meters.
/// @throw         Does not throw (see notes).
///
/// @note If the denominator (den) is zero then funny things could happen; this
///       however should **never** occur for any reference ellipsoid.
///
/// @see "Physical Geodesy", pg. 194
/// @see https://en.wikipedia.org/wiki/Earth_radius
///
template<ellipsoid E>
    double N(double lat) noexcept
{
    return core::N(lat, ellipsoid_traits<E>::a, semi_minor<E>());
}

/// @brief Compute the meridional radii of curvature at a given latitude (on a 
///        reference ellipsoid).
///
/// @param[in] lat The latitude in radians.
/// @tparam    E   The reference ellipsoid (i.e. one of ngpt::ellipsoid).
/// @return        The meridional radius of curvature in meters.
/// @throw         Does not throw (see notes).
///
/// @see https://en.wikipedia.org/wiki/Earth_radius
///
template<ellipsoid E>
    double M(double lat) noexcept
{
    return core::M(lat, ellipsoid_traits<E>::a, semi_minor<E>());
}

class Ellipsoid
{
public:
    Ellipsoid(ellipsoid e) noexcept
    : __ell(e)
    {
        switch (__ell) {
            case ellipsoid::grs80:
                __a = ellipsoid_traits<ellipsoid::grs80>::a;
                __f = ellipsoid_traits<ellipsoid::grs80>::f;
                break;
            case ellipsoid::wgs84:
                __a = ellipsoid_traits<ellipsoid::wgs84>::a;
                __f = ellipsoid_traits<ellipsoid::wgs84>::f;
                break;
            case ellipsoid::pz90:
                __a = ellipsoid_traits<ellipsoid::pz90>::a;
                __f = ellipsoid_traits<ellipsoid::pz90>::f;
                break;
        }
    }

    double
    semi_major() const noexcept
    { return __a; }

    double
    flattening() const noexcept
    { return __f; }

    double
    eccentricity_squared() const noexcept
    { return core::eccentricity_squared(__f); }

    double
    semi_minor() const noexcept
    { return core::semi_minor(__a, __f); }

    double
    N(double lat) const noexcept
    { return core::N(lat, __a, this->semi_minor()); }
    
    double
    M(double lat) const noexcept
    { return core::M(lat, __a, this->semi_minor()); }
private:
    ellipsoid __ell;
    double    __a,
              __f;
}; // class Ellipsoid

} // end namespace ngpt

#endif
