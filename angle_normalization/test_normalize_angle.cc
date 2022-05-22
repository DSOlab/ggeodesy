#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <limits>
#include <random>

using namespace std::chrono;
constexpr const double tolerance = std::numeric_limits<double>::epsilon();

constexpr const double pi = M_PI;
constexpr const double twopi = 2e0 * M_PI;
constexpr const int num_tests = 10000;

std::random_device rd;  // obtain a random number from hardware
std::mt19937 gen(rd()); // seed the generator
std::uniform_real_distribution<double> distr(-10e0 * pi,
                                             10e0 * pi); // define the range

double alg1(double a) noexcept {
  a = std::fmod(a, twopi);
  return (a < 0e0) ? (a + twopi) : a;
}

double alg2(double a) noexcept {
  a = std::fmod(a, twopi);
  const double r[] = {a, a + twopi};
  return r[a < 0e0];
}

double alg3(double a) noexcept {
  return a - twopi * std::floor(a * (1e0 / twopi));
}

double Frac(double x) noexcept { return x - std::floor(x); };
double alg4(double x) noexcept { return twopi * Frac(x / twopi); }

int dummy_op(double a1, double *a2, int sz) noexcept {
  int gtz = 0;
  for (int i = 0; i < sz; i++)
    if (a1 < a2[i])
      ++gtz;
  return gtz;
}

#ifdef ALG1
const char *algorithm_name = "alg1";
#elif ALG2
const char *algorithm_name = "alg2";
#elif ALG3
const char *algorithm_name = "alg3";
#elif ALG4
const char *algorithm_name = "alg4";
#endif

int main() {
  for (int i = 0; i < num_tests; i++) {
    double a = distr(gen);
    double r[4];
    r[0] = alg1(a);
    r[1] = alg2(a);
    r[2] = alg3(a);
    r[3] = alg4(a);
    for (int j = 0; j < 2; j++)
      assert(r[j] == r[j + 1]);
  }
  // printf("All algorithms produce identical results.\n");

  double rand[60];
  for (int i = 0; i < 60; i++)
    rand[i] = distr(gen);

  unsigned long gtz = 0; // dummy var, avoid too much optimization ...
  auto start = high_resolution_clock::now();
  for (int i = 0; i < num_tests; i++) {
    double a = distr(gen);
#ifdef ALG1
    a = alg1(a);
#elif ALG2
    a = alg2(a);
#elif ALG3
    a = alg3(a);
#elif ALG4
    a = alg4(a);
#endif
    for (int k = 0; k < 60; k = k + 2)
      rand[k] = a + distr(gen);
    gtz+=dummy_op(a, rand, 60);
  }
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<microseconds>(stop - start);
  printf(
      "Testing performance of algorithm %s for %d tests; Running time: %lu dummy=%lu\n",
      algorithm_name, num_tests, duration.count(),gtz);

  return 0;
}
