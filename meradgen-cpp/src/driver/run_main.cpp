#include "meradgen_api.hpp"
#include "meradgen_globals.hpp"
#include <cmath>
#include <iostream>
#include <random>

int main() {
  using namespace meradgen;

  const double elab = 10.6;
  const double m = 0.511e-3;
  const double m2 = 0.261112e-6;
  const double pi = std::atan(1.0) * 4.0;

  itest = 0;
  s = 2.0 * (elab * m + m2);
  const double ecm = std::sqrt(2.0 * m * (elab + m)) / 2.0;
  const double pcm = std::sqrt(ecm * ecm - m2);

  merad_init(elab);

  const double thetacm = 90.0;
  const double phi = 10.0;
  const double pl = -1.0;

  float vpgen[4];
  vpgen[3] = 0.0f;  // E (virtual photon)
  vpgen[0] = static_cast<float>(-pcm * std::sin(thetacm * pi / 180.0) * std::cos(phi * pi / 180.0));
  vpgen[1] = static_cast<float>(-pcm * std::sin(thetacm * pi / 180.0) * std::sin(phi * pi / 180.0));
  vpgen[2] = static_cast<float>(pcm * (1.0 - std::cos(thetacm * pi / 180.0)));

  // Number of events to generate
  const int n = 100;

  //////////////////////////////////////////////////////////////////////////////////////////////////////
  // RANDOM NUMBER GENERATION HANDLED HERE IN DRIVER. NO NEED TO BE CONCERNED WITH THIS IN THE MAIN CODE.
  // CHOICE OF SEED AND RANDOM NUMBER GENERATION METHOD IS LEFT TO THE USER. THIS IS FOR EXAMPLE ONLY.

  // Seed for the random number generator
  const uint32_t seed = 12345u;
  // Initialize the random number generator with the seed -- from <random> library
  std::mt19937 rng(seed);
  // Set the random number generator to a uniform distribution between 0 and 1
  std::uniform_real_distribution<double> uni01(0.0, 1.0);
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////

  // SAMPLE RUN -- SAME RANDOMS DIFFERENT POLARIZATION
  for (int i = 0; i < n; i++) {
    std::cout << "EVENT " << i << " START\n";
    double rand4[4] = {uni01(rng), uni01(rng), uni01(rng), uni01(rng)};
    std::cout << "random numbers: " << rand4[0] << " " << rand4[1] << " " << rand4[2] << " " << rand4[3] << "\n";
    meradgen::meradgen(pl, vpgen, rand4);
    std::cout << "polarization: " << pl << "\n";
    std::cout << "vprad: ";
    std::cout << vprad[0] << " " << vprad[1] << " " << vprad[2] << " " << vprad[3] << "\n";
    std::cout << "phrad: ";
    std::cout << phirad[0] << " " << phirad[1] << " " << phirad[2] << " " << phirad[3] << "\n";
    std::cout << "weight: " << weight << "\n";
    meradgen::meradgen(-pl, vpgen, rand4);
    std::cout << "polarization: " << -pl << "\n";
    std::cout << "vprad: ";
    std::cout << vprad[0] << " " << vprad[1] << " " << vprad[2] << " " << vprad[3] << "\n";
    std::cout << "phrad: ";
    std::cout << phirad[0] << " " << phirad[1] << " " << phirad[2] << " " << phirad[3] << "\n";
    std::cout << "weight: " << weight << "\n";
  }
  
  return 0;
}
