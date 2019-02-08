
namespace core
{
    /// @brief Check if a given UTM zone is valid (aka is in range [1,60])
    /// @param[in] zone UTM zone to check
    /// @return false means that the zone is outside range; true otherwise
    /// @see http://www.jaworski.ca/utmzones.htm
    bool
    check_valid_utm_zone(int zone)
    noexcept
    {return zone > 0 && zone <= 60;}
} // namespace core

/// @brief Constants and coefficients to us in utm transformations
namespace utm_constants
{
    constexpr double K0 = 0.9996e0;

    constexpr double E = 0.00669438e0;
    constexpr double E2 = E * E;
    constexpr double E3 = E2 * E;
    constexpr double E_P2 = E / (1e0 - E);

    constexpr double SQRT_E = std::sqrt(1e0 - E);
    constexpr double _E = (1e0 - SQRT_E) / (1e0 + SQRT_E);
    constexpr double _E2 = _E * _E;
    constexpr double _E3 = _E2 * _E;
    constexpr double _E4 = _E3 * _E;
    constexpr double _E5 = _E4 * _E;

    constexpr double M1 = (1e0 - E / 4e0 - 3e0 * E2 / 64e0 - 5e0 * E3 / 256e0);
    constexpr double M2 = (3e0 * E / 8e0 + 3e0 * E2 / 32e0 + 45e0 * E3 / 1024e0);
    constexpr double M3 = (15e0 * E2 / 256e0 + 45e0 * E3 / 1024e0);
    constexpr double M4 = (35e0 * E3 / 3072e0);

    constexpr double P2 = (3e0 / 2e0 * _E - 27e0 / 32e0 * _E3 + 269e0 / 512e0 * _E5);
    constexpr double P3 = (21e0 / 16e0 * _E2 - 55e0 / 32e0 * _E4);
    constexpr double P4 = (151e0 / 96e0 * _E3 - 417e0 / 128e0 * _E5);
    constexpr double P5 = (1097e0 / 512e0 * _E4);
} // namespace utm_constants

/// @brief convert an UTM coordinate into Latitude and Longitude
///
/// @param[in] easting  Easting value of UTM coordinate (meters)
/// @param[in] northing Northing value of UTM coordinate (meters)
/// @param[in] zone_number Zone Number is represented with global map numbers of an UTM Zone
/// @see https://github.com/Turbo87/utm/blob/master/utm/conversion.py
void
utm2ell(double easting, double northing, int zone, double& lon, double& lat, 
    bool northern=false)
{
    using namespace utm_constants;

    if (!core::check_valid_utm_zone(zone)) {
        throw std::out_of_range("[ERROR] Invalid UTM zone");
    }


    double x = easting - 500000e0;
    double y = northing;
    if (!northern) y -= 10000000e0;
    
    double m = y / K0;
    double mu = m / (R * M1);

    double p_rad = (mu +
             P2 * std::sin(2e0 * mu) +
             P3 * std::sin(4e0 * mu) +
             P4 * std::sin(6e0 * mu) +
             P5 * std::sin(8e0 * mu));

    double p_sin  = std::sin(p_rad);
    double p_sin2 = p_sin * p_sin;

    double p_cos = std::cos(p_rad);

    double p_tan  = p_sin / p_cos;
    double p_tan2 = p_tan * p_tan;
    double p_tan4 = p_tan2 * p_tan2;

    double ep_sin = 1e0 - E * p_sin2;
    double ep_sin_sqrt = std::sqrt(1e0 - E * p_sin2);

    double n = R / ep_sin_sqrt;
    double r = (1e0 - E) / ep_sin;

    double c = _E * (p_cos*p_cos);
    double c2 = c * c;

    double d  = x / (n * K0);
    double d2 = d * d;
    double d3 = d2 * d;
    double d4 = d3 * d;
    double d5 = d4 * d;
    double d6 = d5 * d;

    latitude = (p_rad - (p_tan / r) *
        (d2/2e0 -
         d4/24e0 *  (5e0+3e0*p_tan2+10e0*c-4e0*c2-9e0*E_P2)) +
         d6/720e0 * (61e0+90e0*p_tan2+298e0*c+45e0*p_tan4-252e0*E_P2-3e0*c2));

    longitude = (d -
        d3 / 6e0 * (1e0+2e0*p_tan2+c) +
        d5 / 120e0 * (5e0-2e0*c+28e0*p_tan2-3e0*c2+8e0*E_P2+24e0*p_tan4))
        / p_cos;

    return;
}
