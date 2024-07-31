#include "transformations.hpp"
#include <cassert>
#include <cstdio>

using namespace dso;

void touch(CartesianCrdView v) noexcept {
  v.x() -= 1e0;
  v.y() -= 1e0;
  v.z() -= 1e0;
}

int main() {

  Eigen::Matrix<double,3,1> v;
  v(0) = 1e0;
  v(1) = 2e0;
  v(2) = 3e0;

  /* contains its own vector, does not alter v */
  CartesianCrd g1(v);
  assert(g1.x() == 1e0);
  assert(g1.y() == 2e0);
  assert(g1.z() == 3e0);

  /* a view of v; will alter original vector v, but not g1 */
  CartesianCrdView g2(v);
  g2.y() -= 2e0;
  assert(g2.x() == 1e0);
  assert(g2.y() == 0e0);
  assert(g2.z() == 3e0);

  /* alter the instance g1; v and g2 should remain unchanged */
  g1.y() = 4e0;
  assert(g2.x() == 1e0);
  assert(g2.y() == 0e0);
  assert(g2.z() == 3e0);
  assert(v(0) == 1e0);
  assert(v(1) == 0e0);
  assert(v(2) == 3e0);

  /* alter the g2 instance via a function; should not change g1 */
  touch(g2);
  assert(g2.x() == 0e0);
  assert(g2.y() == -1e0);
  assert(g2.z() == 2e0);
  assert(g1.x() == 1e0);
  assert(g1.y() == 4e0);
  assert(g1.z() == 3e0);

  /* alter the g1 instance via a function */
  touch(CartesianCrdView(g1));
  assert(g1.x() == 0e0);
  assert(g1.y() == 3e0);
  assert(g1.z() == 2e0);

  return 0;
}
