// closure_approach_a — Approach A method-validation harness.
//
// Modes:
//   direct  — sample at --pl (default +1); dump lab E / theta + MERADGEN pieces
//   ref0    — sample_reference(P=0); weight_at(+1/-1) when MERADGEN_HAS_WEIGHT_AT
//   paired  — same rand4 at P=+1 and P=0 (bias diagnostic; NOT Approach A)
//
// Prefer raw dumps; use analyze_closure.py for histograms / A.
// See METHOD.md and meradgen-cpp-final-dev/DESIGN_WEIGHT_AT.md.

#include "meradgen_molpol.hpp"
#include "meradgen_globals.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>

namespace {

struct LabParts {
  double scattered[4];
  double recoil[4];
  double photon[4];
};

void boost_cm_to_lab(double elab, double pxpypzE_cm[4]) {
  using meradgen::m;
  using meradgen::s;
  const double gamma = (elab + m) / std::sqrt(s);
  const double beta = std::sqrt(std::max(0.0, 1.0 - 1.0 / (gamma * gamma)));
  const double px = pxpypzE_cm[0];
  const double py = pxpypzE_cm[1];
  const double pz = pxpypzE_cm[2];
  const double e = pxpypzE_cm[3];
  pxpypzE_cm[0] = px;
  pxpypzE_cm[1] = py;
  pxpypzE_cm[2] = gamma * (pz + beta * e);
  pxpypzE_cm[3] = gamma * (e + beta * pz);
}

void reconstruct_lab(const meradgen::MolPolEvent& ev, double elab, LabParts& out) {
  using meradgen::m2;
  using meradgen::s;
  const double ecm = std::sqrt(s) / 2.0;
  const double pcm = std::sqrt(ecm * ecm - m2);

  out.recoil[0] = ev.vprad[0];
  out.recoil[1] = ev.vprad[1];
  out.recoil[2] = ev.vprad[2] - pcm;
  out.recoil[3] = ev.vprad[3] + ecm;

  for (int i = 0; i < 4; ++i)
    out.photon[i] = ev.phirad[i];

  out.scattered[0] = -out.recoil[0] - out.photon[0];
  out.scattered[1] = -out.recoil[1] - out.photon[1];
  out.scattered[2] = -out.recoil[2] - out.photon[2];
  out.scattered[3] = 2.0 * ecm - out.recoil[3] - out.photon[3];

  boost_cm_to_lab(elab, out.scattered);
  boost_cm_to_lab(elab, out.recoil);
  boost_cm_to_lab(elab, out.photon);
}

double lab_theta(const double pxpypzE[4]) {
  const double p = std::sqrt(pxpypzE[0] * pxpypzE[0] + pxpypzE[1] * pxpypzE[1] +
                             pxpypzE[2] * pxpypzE[2]);
  if (!(p > 0.0))
    return 0.0;
  return std::acos(std::max(-1.0, std::min(1.0, pxpypzE[2] / p)));
}

struct Args {
  std::string mode = "paired";
  double elab = 11.0;
  double thetacm_deg = 90.0;
  double phi_deg = 0.0;
  double pl = 1.0;
  double pl_ref = 0.0;
  int n = 500;
  std::uint32_t seed = 1;
  bool radiative_only = false;
  std::string out_path = "kinematics_examination/out/closure_smoke.txt";
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
    if (k == "--mode")
      a.mode = need("--mode");
    else if (k == "--elab")
      a.elab = std::atof(need("--elab"));
    else if (k == "--thetacm")
      a.thetacm_deg = std::atof(need("--thetacm"));
    else if (k == "--phi")
      a.phi_deg = std::atof(need("--phi"));
    else if (k == "--pl")
      a.pl = std::atof(need("--pl"));
    else if (k == "--pl-ref")
      a.pl_ref = std::atof(need("--pl-ref"));
    else if (k == "--n")
      a.n = std::atoi(need("--n"));
    else if (k == "--seed")
      a.seed = static_cast<std::uint32_t>(std::strtoul(need("--seed"), nullptr, 10));
    else if (k == "--out")
      a.out_path = need("--out");
    else if (k == "--radiative")
      a.radiative_only = true;
    else if (k == "-h" || k == "--help") {
      std::cout
          << "closure_approach_a --mode direct|ref0|paired\n"
             "  [--elab GeV] [--thetacm deg] [--phi deg] [--pl P] [--pl-ref P]\n"
             "  [--n N] [--seed S] [--out FILE] [--radiative]\n"
             "\n"
             "  direct : sample at --pl (default +1)\n"
             "  ref0   : sample_reference(--pl-ref); weight_at(+1/-1) if available\n"
             "  paired : same rand4 at P=+1 and P=0 (bias diagnostic)\n"
             "\n"
             "See kinematics_examination/METHOD.md\n";
      std::exit(0);
    } else {
      std::cerr << "unknown arg: " << k << "\n";
      return false;
    }
  }
  if (a.mode != "direct" && a.mode != "ref0" && a.mode != "paired") {
    std::cerr << "--mode must be direct, ref0, or paired\n";
    return false;
  }
  return true;
}

void ensure_parent_dir(const std::string& path) {
  const auto slash = path.find_last_of('/');
  if (slash == std::string::npos)
    return;
  const std::string dir = path.substr(0, slash);
  const std::string cmd = "mkdir -p '" + dir + "'";
  if (std::system(cmd.c_str()) != 0)
    std::cerr << "warning: mkdir failed for " << dir << "\n";
}

void write_lab(std::ostream& out, const char* tag, const LabParts& lab) {
  out << ' ' << tag << "_E_scat " << lab.scattered[3]
      << ' ' << tag << "_E_rec " << lab.recoil[3]
      << ' ' << tag << "_E_gam " << lab.photon[3]
      << ' ' << tag << "_th_scat " << lab_theta(lab.scattered)
      << ' ' << tag << "_th_rec " << lab_theta(lab.recoil);
}

void write_pieces(std::ostream& out, const char* tag, const meradgen::MolPolEvent& e) {
  out << ' ' << tag << "_ich " << e.ich
      << ' ' << tag << "_v " << e.vgen
      << ' ' << tag << "_t1 " << e.t1gen
      << ' ' << tag << "_z " << e.zgen
      << ' ' << tag << "_w " << e.weight
      << ' ' << tag << "_xs0 " << e.xs0
      << ' ' << tag << "_sirad " << e.sirad
      << ' ' << tag << "_sinonr " << e.sinonr;
}

#if defined(MERADGEN_HAS_WEIGHT_AT) && MERADGEN_HAS_WEIGHT_AT
void write_weight_pieces(std::ostream& out, const char* tag,
                         const meradgen::WeightPieces& w) {
  out << ' ' << tag << "_xs0 " << w.xs0
      << ' ' << tag << "_sinonr " << w.sinonr
      << ' ' << tag << "_fsir_z " << w.fsir_z
      << ' ' << tag << "_dens " << w.dens
      << ' ' << tag << "_dens_ref " << w.dens_ref
      << ' ' << tag << "_lr " << w.lr
      << ' ' << tag << "_sitot_ref " << w.sitot_ref
      << ' ' << tag << "_w " << w.weight;
}
#endif

}  // namespace

int main(int argc, char** argv) {
  Args args;
  if (!parse_args(argc, argv, args))
    return 2;
  if (args.n <= 0) {
    std::cerr << "--n must be positive\n";
    return 2;
  }

  ensure_parent_dir(args.out_path);
  std::ofstream out(args.out_path);
  if (!out) {
    std::cerr << "cannot write " << args.out_path << "\n";
    return 1;
  }

  const bool has_weight_at =
#if defined(MERADGEN_HAS_WEIGHT_AT) && MERADGEN_HAS_WEIGHT_AT
      true;
#else
      false;
#endif

  using meradgen::pi;
  const double thetacm = args.thetacm_deg * pi / 180.0;
  const double phi = args.phi_deg * pi / 180.0;

  meradgen::itest = 0;
  meradgen::merad_init(args.elab);
  double vp[4];
  meradgen::vpgen_from_angles(args.elab, thetacm, phi, vp);

  out << "# closure_approach_a dump\n";
  out << "# mode " << args.mode
      << " elab_GeV " << args.elab
      << " thetacm_deg " << args.thetacm_deg
      << " phi_deg " << args.phi_deg;
  if (args.mode == "direct")
    out << " pl_sample " << args.pl;
  if (args.mode == "ref0")
    out << " pl_ref " << args.pl_ref;
  out << "\n";
  out << "# target_events " << args.n << " seed " << args.seed;
  if (args.radiative_only)
    out << " radiative_only 1";
  out << " has_weight_at " << (has_weight_at ? 1 : 0) << "\n";
  out << "# vpgen_pxpypzE " << std::setprecision(16) << vp[0] << ' ' << vp[1]
      << ' ' << vp[2] << ' ' << vp[3] << "\n";
  out << "# lab (px,py,pz,E) GeV; th_* radians; WeightPieces per DESIGN_WEIGHT_AT.md\n";
  out << "# one record per accepted event (space-separated key/value tokens)\n";
  if (!has_weight_at && args.mode == "ref0") {
    out << "# NOTE: weight_at not compiled in — ref0 dumps pl_ref sample only; "
           "see METHOD.md\n";
  }

  std::mt19937 rng(args.seed);
  std::uniform_real_distribution<double> U(0.0, 1.0);

  std::size_t written = 0;
  std::size_t tried = 0;
  const std::size_t max_tries =
      args.radiative_only ? static_cast<std::size_t>(args.n) * 1000u
                          : static_cast<std::size_t>(args.n) * 20u;

  while (static_cast<int>(written) < args.n && tried < max_tries) {
    const double rand4[4] = {U(rng), U(rng), U(rng), U(rng)};
    ++tried;

    if (args.mode == "direct") {
      meradgen::MolPolEvent e{};
      const bool ok = meradgen::generate(args.pl, vp, rand4, e);
      if (!ok)
        continue;
      if (args.radiative_only && e.ich != 1)
        continue;
      LabParts lab{};
      reconstruct_lab(e, args.elab, lab);
      ++written;
      out << "event " << written << " r0 " << std::scientific << std::setprecision(16)
          << rand4[0] << " r1 " << rand4[1] << " r2 " << rand4[2] << " r3 " << rand4[3]
          << " ok 1";
      write_pieces(out, "S", e);
      write_lab(out, "S", lab);
      out << '\n';
    } else if (args.mode == "ref0") {
      meradgen::MolPolEvent e0{};
#if defined(MERADGEN_HAS_WEIGHT_AT) && MERADGEN_HAS_WEIGHT_AT
      const bool ok = meradgen::sample_reference(vp, rand4, e0, args.pl_ref);
#else
      const bool ok = meradgen::generate(args.pl_ref, vp, rand4, e0);
#endif
      if (!ok)
        continue;
      if (args.radiative_only && e0.ich != 1)
        continue;
      LabParts lab{};
      reconstruct_lab(e0, args.elab, lab);
      ++written;
      out << "event " << written << " r0 " << std::scientific << std::setprecision(16)
          << rand4[0] << " r1 " << rand4[1] << " r2 " << rand4[2] << " r3 " << rand4[3]
          << " ok 1";
      write_pieces(out, "P0", e0);
      write_lab(out, "P0", lab);
#if defined(MERADGEN_HAS_WEIGHT_AT) && MERADGEN_HAS_WEIGHT_AT
      meradgen::WeightPieces wp{}, wm{};
      const bool okp = meradgen::weight_at(+1.0, vp, e0, args.pl_ref, wp);
      const bool okm = meradgen::weight_at(-1.0, vp, e0, args.pl_ref, wm);
      out << " wa_ok_pp " << (okp ? 1 : 0) << " wa_ok_pm " << (okm ? 1 : 0);
      if (okp)
        write_weight_pieces(out, "Wp", wp);
      if (okm)
        write_weight_pieces(out, "Wm", wm);
#endif
      out << '\n';
    } else {  // paired
      meradgen::MolPolEvent ep{}, e0{};
      const bool okp = meradgen::generate(+1.0, vp, rand4, ep);
      const bool ok0 = meradgen::generate(0.0, vp, rand4, e0);
      if (!okp || !ok0)
        continue;
      if (args.radiative_only && !(ep.ich == 1 && e0.ich == 1))
        continue;
      LabParts labp{}, lab0{};
      reconstruct_lab(ep, args.elab, labp);
      reconstruct_lab(e0, args.elab, lab0);
      ++written;
      out << "event " << written << " r0 " << std::scientific << std::setprecision(16)
          << rand4[0] << " r1 " << rand4[1] << " r2 " << rand4[2] << " r3 " << rand4[3]
          << " ok_pp 1 ok_p0 1";
      write_pieces(out, "Pp", ep);
      write_lab(out, "Pp", labp);
      write_pieces(out, "P0", e0);
      write_lab(out, "P0", lab0);
      out << " dE_scat " << (labp.scattered[3] - lab0.scattered[3])
          << " dE_rec " << (labp.recoil[3] - lab0.recoil[3]) << '\n';
    }
  }

  if (static_cast<int>(written) < args.n) {
    std::cerr << "only wrote " << written << " of " << args.n << " (tried " << tried
              << " quads)\n";
    return 1;
  }

  std::cerr << "wrote " << args.out_path << " (" << written << " events, tried "
            << tried << " quads, has_weight_at=" << (has_weight_at ? 1 : 0) << ")\n";
  return 0;
}
