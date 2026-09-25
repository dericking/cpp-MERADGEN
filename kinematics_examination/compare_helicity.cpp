// compare_helicity — same rand4 at P=+1, 0, -1; one dump file.
// Default: elab=11 GeV, thetacm=90°, phi=0° (CM).
//
// usage:
//   compare_helicity [--elab GeV] [--thetacm deg] [--phi deg]
//                    [--n N] [--seed S] [--quads FILE] [--out FILE]
//                    [--radiative]

#include "meradgen_molpol.hpp"
#include "meradgen_globals.hpp"

#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

namespace {

struct Args {
  double elab = 11.0;
  double thetacm_deg = 90.0;
  double phi_deg = 0.0;
  int n = 200;
  std::uint32_t seed = 1;
  bool radiative_only = false;
  std::string quads_path;
  std::string out_path = "kinematics_examination/out/pol_compare_11GeV.txt";
};

bool parse_args(int argc, char** argv, Args& a) {
  for (int i = 1; i < argc; ++i) {
    const std::string k = argv[i];
    auto need = [&](const char* name) -> const char* {
      if (i + 1 >= argc) {
        std::cerr << "missing value for " << name << "\n";
        std::exit(2);
      }
      return argv[++i];
    };
    if (k == "--elab")
      a.elab = std::atof(need("--elab"));
    else if (k == "--thetacm")
      a.thetacm_deg = std::atof(need("--thetacm"));
    else if (k == "--phi")
      a.phi_deg = std::atof(need("--phi"));
    else if (k == "--n")
      a.n = std::atoi(need("--n"));
    else if (k == "--seed")
      a.seed = static_cast<std::uint32_t>(std::strtoul(need("--seed"), nullptr, 10));
    else if (k == "--quads")
      a.quads_path = need("--quads");
    else if (k == "--out")
      a.out_path = need("--out");
    else if (k == "--radiative")
      a.radiative_only = true;
    else if (k == "-h" || k == "--help") {
      std::cout << "compare_helicity [--elab GeV] [--thetacm deg] [--phi deg]\n"
                   "  [--n N] [--seed S] [--quads FILE] [--out FILE] [--radiative]\n"
                   "  --radiative: keep only events with ich==1 for Pp, P0, and Pm\n"
                   "               (with --n, draw until N such events are written)\n";
      std::exit(0);
    } else {
      std::cerr << "unknown arg: " << k << "\n";
      return false;
    }
  }
  return true;
}

std::vector<std::array<double, 4>> load_quads_file(const std::string& path) {
  std::vector<std::array<double, 4>> q;
  std::ifstream in(path);
  if (!in) {
    std::cerr << "cannot open quads: " << path << "\n";
    std::exit(1);
  }
  std::array<double, 4> row{};
  while (in >> row[0] >> row[1] >> row[2] >> row[3])
    q.push_back(row);
  return q;
}

void write_event_block(std::ostream& out, const char* tag,
                       const meradgen::MolPolEvent& e) {
  out << tag << "_ich " << e.ich
      << " " << tag << "_v " << std::scientific << std::setprecision(16) << e.vgen
      << " " << tag << "_t1 " << e.t1gen
      << " " << tag << "_z " << e.zgen
      << " " << tag << "_w " << e.weight
      << " " << tag << "_xs0 " << e.xs0
      << " " << tag << "_sirad " << e.sirad
      << " " << tag << "_sinonr " << e.sinonr;
  out << " " << tag << "_vprad";
  for (int i = 0; i < 4; ++i)
    out << ' ' << e.vprad[i];
  out << " " << tag << "_phirad";
  for (int i = 0; i < 4; ++i)
    out << ' ' << e.phirad[i];
}

bool all_radiative(const meradgen::MolPolEvent& a, const meradgen::MolPolEvent& b,
                   const meradgen::MolPolEvent& c) {
  return a.ich == 1 && b.ich == 1 && c.ich == 1;
}

} // namespace

int main(int argc, char** argv) {
  Args args;
  if (!parse_args(argc, argv, args))
    return 2;

  if (args.n <= 0) {
    std::cerr << "--n must be positive\n";
    return 2;
  }

  // Ensure output directory exists (best-effort; caller may mkdir).
  {
    const auto slash = args.out_path.find_last_of('/');
    if (slash != std::string::npos) {
      const std::string dir = args.out_path.substr(0, slash);
      const std::string cmd = "mkdir -p '" + dir + "'";
      if (std::system(cmd.c_str()) != 0) {
        std::cerr << "warning: mkdir failed for " << dir << "\n";
      }
    }
  }

  std::ofstream out(args.out_path);
  if (!out) {
    std::cerr << "cannot write " << args.out_path << "\n";
    return 1;
  }

  using meradgen::pi;
  const double thetacm = args.thetacm_deg * pi / 180.0;
  const double phi = args.phi_deg * pi / 180.0;

  meradgen::merad_init(args.elab);
  double vp[4];
  meradgen::vpgen_from_angles(args.elab, thetacm, phi, vp);

  out << "# compare_helicity dump\n";
  out << "# elab_GeV " << args.elab
      << " thetacm_deg " << args.thetacm_deg
      << " phi_deg " << args.phi_deg << "\n";
  out << "# target_events " << args.n;
  if (args.radiative_only)
    out << " radiative_only 1";
  if (args.quads_path.empty())
    out << " seed " << args.seed;
  else
    out << " quads_file " << args.quads_path;
  out << "\n";
  out << "# vpgen_pxpypzE " << std::setprecision(16) << vp[0] << ' ' << vp[1] << ' '
      << vp[2] << ' ' << vp[3] << "\n";
  out << "# P = pb*pt; blocks: Pp (P=+1), P0 (P=0), Pm (P=-1); same rand4\n";
  out << "# one record per quad (space-separated key/value tokens)\n";

  auto emit = [&](std::size_t event_id, const double rand4[4],
                  const meradgen::MolPolEvent& ep, const meradgen::MolPolEvent& e0,
                  const meradgen::MolPolEvent& em, bool okp, bool ok0, bool okm) {
    out << "event " << event_id
        << " r0 " << std::setprecision(16) << rand4[0]
        << " r1 " << rand4[1]
        << " r2 " << rand4[2]
        << " r3 " << rand4[3]
        << " ok_pp " << (okp ? 1 : 0)
        << " ok_p0 " << (ok0 ? 1 : 0)
        << " ok_pm " << (okm ? 1 : 0)
        << ' ';
    write_event_block(out, "Pp", ep);
    out << ' ';
    write_event_block(out, "P0", e0);
    out << ' ';
    write_event_block(out, "Pm", em);
    out << '\n';
  };

  std::size_t written = 0;
  std::size_t tried = 0;

  if (!args.quads_path.empty()) {
    const auto quads = load_quads_file(args.quads_path);
    for (std::size_t iq = 0; iq < quads.size() && static_cast<int>(written) < args.n; ++iq) {
      const double rand4[4] = {quads[iq][0], quads[iq][1], quads[iq][2], quads[iq][3]};
      meradgen::MolPolEvent ep{}, e0{}, em{};
      const bool okp = meradgen::generate(+1.0, vp, rand4, ep);
      const bool ok0 = meradgen::generate(0.0, vp, rand4, e0);
      const bool okm = meradgen::generate(-1.0, vp, rand4, em);
      ++tried;
      if (args.radiative_only && !all_radiative(ep, e0, em))
        continue;
      ++written;
      emit(written, rand4, ep, e0, em, okp, ok0, okm);
    }
  } else {
    std::mt19937 rng(args.seed);
    std::uniform_real_distribution<double> U(0.0, 1.0);
    const std::size_t max_tries =
        args.radiative_only ? static_cast<std::size_t>(args.n) * 1000u
                            : static_cast<std::size_t>(args.n);
    while (static_cast<int>(written) < args.n && tried < max_tries) {
      const double rand4[4] = {U(rng), U(rng), U(rng), U(rng)};
      meradgen::MolPolEvent ep{}, e0{}, em{};
      const bool okp = meradgen::generate(+1.0, vp, rand4, ep);
      const bool ok0 = meradgen::generate(0.0, vp, rand4, e0);
      const bool okm = meradgen::generate(-1.0, vp, rand4, em);
      ++tried;
      if (args.radiative_only && !all_radiative(ep, e0, em))
        continue;
      ++written;
      emit(written, rand4, ep, e0, em, okp, ok0, okm);
    }
  }

  if (static_cast<int>(written) < args.n) {
    std::cerr << "only wrote " << written << " of " << args.n
              << " (tried " << tried << " quads)\n";
    return 1;
  }

  std::cerr << "wrote " << args.out_path << " (" << written << " events, tried "
            << tried << " quads)\n";
  return 0;
}
