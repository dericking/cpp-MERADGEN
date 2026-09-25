#include "meradgen_globals.hpp"
#include <cmath>

namespace meradgen {

    /////////////////////////////////////////////////////////////////////////////////////
    // Fixed constants. Live values are PDG 2024/25 (production tree).
    // FORTRAN literals stay commented; they remain live in meradgen-cpp/.
    /////////////////////////////////////////////////////////////////////////////////////

    const double pi     = std::acos(-1.0);          // [source C++ idiom]
    // const double pi  = std::atan(1.0) * 4.0;     // [source MERADGEN Fortran merad_init]

    const double alfa   = 0.00729735256;            // fine structure constant [source PDG 2024/25]
    // const double alfa = 0.729735e-2;             // [source MERADGEN Fortran]

    const double m      = 0.51099895e-3;            // electron mass (GeV) [source PDG 2024/25]
    // const double m   = 0.511000e-3;              // [source MERADGEN Fortran]

    const double m2     = m * m;                    // electron mass squared [source PDG, derived]
    // const double m2  = 0.261112e-6;              // [source MERADGEN Fortran; not equal to m*m]

    const double barn   = 0.3893793721e6;             // (hbar*c)^2 (GeV^2 ub) [source PDG]
    // const double barn = 0.389379e6;              // [source MERADGEN Fortran]

    // Muon/tau: PDG 2024/25. Live in vacpol am2 = {m2, mu2, tau2}.
    const double mu   = 105.65837e-3;               // muon mass (GeV) [source PDG 2024/25]
    const double mu2  = mu * mu;
    const double tau  = 1776.93e-3;                 // tau mass (GeV) [source PDG 2024/25]
    const double tau2 = tau * tau;

    /////////////////////////////////////////////////////////////////////////////////////
    // Event specific values
    double En = 0., s = 0., als = 0., coeb = 0., coer = 0., Egmin = 0., t=0., pl=0.;
    double vprad[4] = {0}, phirad[4] = {0}, weight = 0., xs0_save = 0.;
    double sirad_out = 0., sinonr_out = 0.;

    int ich = 0, itest = 0; // ich: indicator of radiative event; itest: holds test type

    /////////////////////////////////////////////////////////////////////////////////////
    // Grid variables
    double grv[NV + 1] = {0};
    double grt1[NT1 + 1] = {0};
    double grz[NZ + 1] = {0};
    
    double vmin  = 0., vmax  = 0., vgen  = 0.;
    double t1min = 0., t1max = 0., t1gen = 0.;
    double zmin  = 0., zmax  = 0., zgen  = 0.;

    double az  = 0., bz  = 0., cz  = 0.;
    double az1 = 0., bz1 = 0., cz1 = 0.;
    double az2 = 0., bz2 = 0., cz2 = 0.;

    double bin[NBIN + 1] = {0};
    double argbin[NBIN + 1] = {0};
    double sigbin[NBIN + 1] = {0};
    double step = 0.;

    /////////////////////////////////////////////////////////////////////////////////////
    // dcanc() variables
    double xxs = 0., xxt = 0., xxu = 0.;

}
