#include "transformations.hpp"

using namespace dso;

int main() {

  Eigen::Matrix<double,3,1> v;
  GeodeticCrd g1;
  GeodeticCrdView g2(v);
  GeodeticCrdConstView g3(v);

  return 0;
}
