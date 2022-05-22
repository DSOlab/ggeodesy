# Normalizing angles in range [0, 2π]

Obviously, we can perform this operation in a series of different 
implementation approaches. Here are 4:

```
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

double frac(double x) noexcept { return x - std::floor(x); };
double alg4(double x) noexcept { return twopi * frac(x / twopi); }
```

Note that the simple `std::fmod(angle, twopi)` will not yield correct results
for negative angles.

This folder contains a C++ source code file to test the above implementations 
and two bash scripts, one used for building the source and the other for 
actually running/performing the tests.

It seems that no significant performance gains can be obtained by adopting any 
of the above implementations when using a no-compiler-optimization approach 
(aka via `-O0`). When using an optimization level `-O2`, with the GNU compiler, 
the second implementation (`alg2`, the branchless version) seems to be a tiny 
bit faster.

Note however that the results are heavily dependent on the compiler; here i 
test using the GNU compiler (g++) and the clang (clang++).
