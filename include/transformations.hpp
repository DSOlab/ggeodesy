/** @file
 * A list of frequently used geodetic functions.
 */

#ifndef __DSO_GEO_COORDINATE_TRANSFORMATIONS_HPP__
#define __DSO_GEO_COORDINATE_TRANSFORMATIONS_HPP__

#include "core/crd_transformations.hpp"
#include "core/crdtype_warppers.hpp"

namespace dso {

template<typename C=CartesianCrd>
inline SphericalCrd cartesian2spherical(const /*CartesianCrd*/C &v) noexcept {
  static_assert(dso::CoordinateTypeTraits<C>::isCartesian);
  SphericalCrd s;
  cartesian2spherical(v.x(), v.y(), v.z(), s.r(), s.lat(), s.lon());
  return s;
}

template<typename S=SphericalCrd>
inline CartesianCrd spherical2cartesian(const /*SphericalCrd*/S &v) noexcept {
  static_assert(dso::CoordinateTypeTraits<S>::isSpherical);
  CartesianCrd s;
  spherical2cartesian(v.r(), v.lat(), v.lon(), s.x(), s.y(), s.z());
  return s;
}

template <ellipsoid E, typename G=GeodeticCrd>
CartesianCrd geodetic2cartesian(const /*GeodeticCrd*/G &v) noexcept {
  static_assert(dso::CoordinateTypeTraits<G>::isGeodetic);
  CartesianCrd s;
  geodetic2cartesian<E>(v.lat(), v.lon(), v.hgt(), s.x(), s.y(), s.z());
  return s;
}

template <ellipsoid E, typename C=CartesianCrd>
GeodeticCrd cartesian2geodetic(const /*CartesianCrd*/C &v) noexcept {
  static_assert(dso::CoordinateTypeTraits<C>::isCartesian);
  GeodeticCrd s;
  cartesian2geodetic<E>(v.x(), v.y(), v.z(), s.lat(), s.lon(), s.hgt());
  return s;
}
} /* namespace dso */

#endif
