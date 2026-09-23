#include "meradgen_api.hpp"
#include "meradgen_globals.hpp"

#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

static void write_f010(std::ostream& out, double x) {
  char buf[32];
  std::snprintf(buf, sizeof(buf), "%.10f", x);
  std::string s(buf);
  if (s.size() >= 2 && s[0] == '0' && s[1] == '.') {
    s.erase(0, 1);
  }
  out << s;
}

static void write_f32_hex(std::ostream& out, float x) {
  std::uint32_t u = 0;
  std::memcpy(&u, &x, sizeof u);
  out << "0x" << std::hex << std::uppercase << std::setw(8) << std::setfill('0') << u
      << std::dec << std::nouppercase << std::setfill(' ');
}

static void write_vec4(std::ostream& out, const char* tag, const double v[4],
                       const std::string& prec) {
  out << "  " << tag << "{";
  for (int i = 0; i < 4; ++i) {
    if (i) out << ' ';
    const float f = static_cast<float>(v[i]);
    if (prec == "hex") {
      write_f32_hex(out, f);
    } else if (prec == "full") {
      out << std::scientific << std::uppercase << std::setprecision(16) << static_cast<double>(f);
    } else {
      out << std::scientific << std::uppercase << std::setprecision(6) << std::setw(14)
          << static_cast<double>(f);
    }
  }
  out << "}\n";
}

static void write_kin(std::ostream& out, const std::string& prec) {
  using namespace meradgen;
  out << "  KIN{" << ich << ' ';
  if (prec == "es14") {
    out << std::scientific << std::uppercase << std::setprecision(6) << vgen << ' ' << t1gen << ' '
        << zgen;
  } else {
    out << std::scientific << std::uppercase << std::setprecision(16) << vgen << ' ' << t1gen << ' '
        << zgen;
  }
  out << "}\n";
}

int main(int argc, char** argv) {
  using namespace meradgen;
  if (argc < 3) {
    std::cerr << "usage: driver_parity quads.txt output.txt [es14|full|hex] "
                 "[max_rad] [max_calls] [elab thetacm phi pl]\n";
    return 1;
  }
  const std::string prec = (argc >= 4) ? argv[3] : "es14";
  if (prec != "es14" && prec != "full" && prec != "hex") {
    std::cerr << "unknown prec: " << prec << "\n";
    return 1;
  }
  const int max_rad = (argc >= 5) ? std::atoi(argv[4]) : 0;
  const int max_calls = (argc >= 6) ? std::atoi(argv[5]) : 0;

  float elab_f = 45.f;
  double thetacm = 90.0;
  double phi = 10.0;
  double pl = -1.0;
  if (argc >= 7) {
    elab_f = static_cast<float>(std::atof(argv[6]));
  }
  if (argc >= 8) {
    thetacm = std::atof(argv[7]);
  }
  if (argc >= 9) {
    phi = std::atof(argv[8]);
  }
  if (argc >= 10) {
    pl = std::atof(argv[9]);
  }
  const double elab = static_cast<double>(elab_f);
  itest = 0;
  merad_init(elab);

  const double ecm = std::sqrt(2.0 * m * (elab + m)) / 2.0;
  const double pcm = std::sqrt(ecm * ecm - m2);

  float vpgen[4];
  vpgen[3] = 0.0f;
  vpgen[0] = static_cast<float>(-pcm * std::sin(thetacm * pi / 180.0) *
                                std::cos(phi * pi / 180.0));
  vpgen[1] = static_cast<float>(-pcm * std::sin(thetacm * pi / 180.0) *
                                std::sin(phi * pi / 180.0));
  vpgen[2] = static_cast<float>(pcm * (1.0 - std::cos(thetacm * pi / 180.0)));

  std::ifstream qin(argv[1]);
  if (!qin) {
    std::cerr << "cannot open quads: " << argv[1] << "\n";
    return 1;
  }
  std::ofstream out(argv[2]);
  if (!out) {
    std::cerr << "cannot open output: " << argv[2] << "\n";
    return 1;
  }

    out << "HEADER START\n";
  out << std::scientific << std::uppercase << std::setprecision(16);
  out << "elab=" << elab << "\n";
  out << "thetacm=" << thetacm << "\n";
  out << "phi=" << phi << "\n";
  out << "pl=" << pl << "\n";
  out << "m=" << m << "\n";
  out << "m2=" << m2 << "\n";
  out << "ecm=" << ecm << "\n";
  out << "pcm=" << pcm << "\n";
  out << "vpgen=" << static_cast<double>(vpgen[0]) << ' ' << static_cast<double>(vpgen[1])
      << ' ' << static_cast<double>(vpgen[2]) << ' ' << static_cast<double>(vpgen[3]) << "\n";
  out << std::defaultfloat;
  out << "mode=" << prec << "\n";
  out << "max_rad=" << max_rad << "\n";
  out << "max_calls=" << max_calls << "\n";
  out << "dump_policy=" << (max_rad > 0 ? "radiative_only" : "all") << "\n";
  out << "HEADER END\n";

  int iev = 0;
  int nrad = 0;
  double r1, r2, r3, r4;
  while (qin >> r1 >> r2 >> r3 >> r4) {
    ++iev;
    const double rand4[4] = {r1, r2, r3, r4};
    meradgen::meradgen(pl, vpgen, rand4);
    if (ich == 1) {
      ++nrad;
    }
    if (max_rad == 0 || ich == 1) {
      out << "EVENT " << iev << " START\n";
      out << "  RANDOM{";
      write_f010(out, r1);
      out << ' ';
      write_f010(out, r2);
      out << ' ';
      write_f010(out, r3);
      out << ' ';
      write_f010(out, r4);
      out << "}\n";
      write_vec4(out, "VPRAD", vprad, prec);
      write_vec4(out, "PHIRAD", phirad, prec);
      if (prec != "es14") {
        write_kin(out, prec);
      }
    }
    if (max_rad > 0 && nrad >= max_rad) {
      break;
    }
    if (max_calls > 0 && iev >= max_calls) {
      break;
    }
  }
  out << std::defaultfloat;
  out << "FOOTER START\n";
  out << "calls=" << iev << "\n";
  out << "radiative=" << nrad << "\n";
  out << "FOOTER END\n";
  std::cerr << "calls=" << iev << " radiative=" << nrad << "\n";
  std::cerr << std::scientific << std::setprecision(8)
            << "kinematics elab,thetacm,phi,pl=" << elab << ' ' << thetacm << ' '
            << phi << ' ' << pl << "\n";
  return 0;
}
