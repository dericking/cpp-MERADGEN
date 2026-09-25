#include "meradgen_api.hpp"
#include "meradgen_globals.hpp"

#include <cmath>
#include <cstdint>
#include <cstring>
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

static void write_vec4(std::ostream& out, const char* tag, const double v[4]) {
  out << "  " << tag << "{";
  for (int i = 0; i < 4; ++i) {
    if (i) out << ' ';
    out << std::scientific << std::uppercase << std::setprecision(16) << v[i];
  }
  out << "}\n";
}

int main(int argc, char** argv) {
  using namespace meradgen;
  if (argc < 3) {
    std::cerr << "usage: meradgen_quads quads.txt output.txt\n";
    return 1;
  }

  const double elab = 45.0;
  itest = 0;
  merad_init(elab);

  const double ecm = std::sqrt(2.0 * m * (elab + m)) / 2.0;
  const double pcm = std::sqrt(ecm * ecm - m2);
  const double thetacm = 90.0;
  const double phi = 10.0;
  const double pl = -1.0;

  double vpgen[4];
  vpgen[3] = 0.0;
  vpgen[0] = -pcm * std::sin(thetacm * pi / 180.0) * std::cos(phi * pi / 180.0);
  vpgen[1] = -pcm * std::sin(thetacm * pi / 180.0) * std::sin(phi * pi / 180.0);
  vpgen[2] = pcm * (1.0 - std::cos(thetacm * pi / 180.0));

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

  int iev = 0;
  double r1, r2, r3, r4;
  while (qin >> r1 >> r2 >> r3 >> r4) {
    ++iev;
    const double rand4[4] = {r1, r2, r3, r4};
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
    meradgen::meradgen(pl, vpgen, rand4);
    write_vec4(out, "VPRAD", vprad);
    write_vec4(out, "PHIRAD", phirad);
    out << "  KIN{" << ich << ' ' << std::scientific << std::uppercase
        << std::setprecision(16) << vgen << ' ' << t1gen << ' ' << zgen << "}\n";
  }
  return 0;
}
