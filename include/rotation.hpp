#ifndef __DSO_GEODESY_ROTATIONS_HPP__
#define __DSO_GEODESY_ROTATIONS_HPP__

#include <cmath>
#include "eigen3/Eigen/Eigen"

namespace dso {
  /** @brief Compare 3x3 rotation matrices.
   *
   * Let P and Q be orthogonal matrices representing two rotations in the same 
   * basis. Let Q^T denote the matrix transpose. The difference rotation 
   * matrix that represents the difference rotation is defined as:
   * R = P * Q^T
   * The distance between rotations represented by rotation matrices P and Q 
   * is the angle of the difference rotation represented by the rotation 
   * matrix R.
   * We can retrieve the angle of the difference rotation from the trace of R
   * trR = 1 + 2 * cosθ
   * or
   * θ = acos( (trR - 1) / 2 )
   *
   * @param[in] P Rotation matrix of size 3x3
   * @param[in] Q Rotation matrix of size 3x3
   * @return    The angle theta in [rad] defined as:
   *            θ = acos( (trR - 1) / 2 ), with R =  P * Q^T
   */
inline double rotation_distance(const Eigen::Matrix<double, 3, 3> &P,
                                const Eigen::Matrix<double, 3, 3> &Q) noexcept {
  return std::acos((P * Q.transpose()).trace() / 2e0);
}
} /* namespace dso */

#endif
