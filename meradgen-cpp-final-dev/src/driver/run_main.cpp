#include "meradgen_api.hpp"
#include "meradgen_globals.hpp"
#include "meradgen_molpol.hpp"
#include <cmath>
#include <iostream>
#include <random>

int main() {
  using namespace meradgen;

  const double elab = 10.6;
  itest = 0;
  merad_init(elab);

  const double thetacm = 90.0 * pi / 180.0;
  const double phi = 10.0 * pi / 180.0;
  const double pl = -1.0;

  double vpgen[4];
  vpgen_from_angles(elab, thetacm, phi, vpgen);

  const int n = 100;
  const uint32_t seed = 12345u;
  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> uni01(0.0, 1.0);

  int n_reject = 0;
  for (int i = 0; i < n; i++) {
    double rand4[4] = {uni01(rng), uni01(rng), uni01(rng), uni01(rng)};
    MolPolEvent plus, minus;
    if (!generate_pair(pl, vpgen, rand4, plus, minus)) {
      n_reject++;
      continue;
    }
    std::cout << "EVENT " << i << " START\n";
    std::cout << "random numbers: " << rand4[0] << " " << rand4[1] << " " << rand4[2] << " "
              << rand4[3] << "\n";
    std::cout << "polarization: " << pl << " ich: " << plus.ich << "\n";
    std::cout << "vprad: " << plus.vprad[0] << " " << plus.vprad[1] << " " << plus.vprad[2] << " "
              << plus.vprad[3] << "\n";
    std::cout << "phrad: " << plus.phirad[0] << " " << plus.phirad[1] << " " << plus.phirad[2]
              << " " << plus.phirad[3] << "\n";
    std::cout << "weight: " << plus.weight << "\n";
    std::cout << "polarization: " << -pl << " ich: " << minus.ich << "\n";
    std::cout << "vprad: " << minus.vprad[0] << " " << minus.vprad[1] << " " << minus.vprad[2]
              << " " << minus.vprad[3] << "\n";
    std::cout << "phrad: " << minus.phirad[0] << " " << minus.phirad[1] << " " << minus.phirad[2]
              << " " << minus.phirad[3] << "\n";
    std::cout << "weight: " << minus.weight << "\n";
  }
  std::cout << "rejected_nonfinite=" << n_reject << "\n";
  return 0;
}
