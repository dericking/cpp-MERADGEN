#include <iostream>
#include <cmath>
#include <vector>
#include <random>

void run(double pol, double enlab, std::vector<double> &virtphoton, 
         std::vector<double> &virtphotonradplus, std::vector<double> &realphotonplus, 
         double &vvarplus, double &zvarplus, double &t1varplus, 
         double &sigmabornplus, double &sigmaradplus, double &sigmanoradplus, 
         std::vector<double> &virtphotonradminus, std::vector<double> &realphotonminus, 
         double &vvarminus, double &zvarminus, double &t1varminus, 
         double &sigmabornminus, double &sigmaradminus, double &sigmanoradminus) {
    
    const double pi = atan(1.0) * 4.0;
    const double m = 0.511000e-3;
    const double m2 = 0.261112e-6;
    double elab = enlab;
    int itest = 0;
    double s = 2.0 * (elab * m + m2);
    double pl = 1.0;
    std::vector<double> vpgen(4);
    std::vector<double> r(4);
    
    for (int it = 0; it < 3; ++it) {
        vpgen[it] = virtphoton[it + 1];
    }
    vpgen[3] = virtphoton[0];

    // Pre generation of random numbers to use them in both p=1 and p=-1
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);
    for (auto &val : r) {
        val = dis(gen);
    }

    // Call to merad_init(elab); // Function not defined in the provided code
    // Call to meradgen(pl, vpgen, r); // Function not defined in the provided code

    // Output (output.inc)
    // vprad = p2 - p1 in CM system
    // phirad - 4-momentum of the real photon
    // Call to meradgen(pl, vpgen, r); // Function not defined in the provided code

    for (int it = 0; it < 3; ++it) {
        virtphotonradplus[it + 1] = vprad[it]; // vprad not defined in the provided code
        realphotonplus[it + 1] = phirad[it];   // phirad not defined in the provided code
    }
    virtphotonradplus[0] = vprad[3]; // vprad not defined in the provided code
    realphotonplus[0] = phirad[3];   // phirad not defined in the provided code
    vvarplus = vgen;                  // vgen not defined in the provided code
    zvarplus = zgen;                  // zgen not defined in the provided code
    t1varplus = t1gen;                // t1gen not defined in the provided code
    sigmabornplus = (sigmar + sinonr) / weight; // sigmar, sinonr, weight not defined in the provided code
    sigmaradplus = sigmar;            // sigmar not defined in the provided code
    sigmanoradplus = sinonr;          // sinonr not defined in the provided code

    // Call to meradgen(-pl, vpgen, r); // Function not defined in the provided code

    for (int it = 0; it < 3; ++it) {
        virtphotonradminus[it + 1] = vprad[it]; // vprad not defined in the provided code
        realphotonminus[it + 1] = phirad[it];   // phirad not defined in the provided code
    }
    virtphotonradminus[0] = vprad[3]; // vprad not defined in the provided code
    realphotonminus[0] = phirad[3];   // phirad not defined in the provided code
    vvarminus = vgen;                  // vgen not defined in the provided code
    zvarminus = zgen;                  // zgen not defined in the provided code
    t1varminus = t1gen;                // t1gen not defined in the provided code
    sigmabornminus = (sigmar + sinonr) / weight; // sigmar, sinonr, weight not defined in the provided code
    sigmaradminus = sigmar;            // sigmar not defined in the provided code
    sigmanoradminus = sinonr;          // sinonr not defined in the provided code
}

