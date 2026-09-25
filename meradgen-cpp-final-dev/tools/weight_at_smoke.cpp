// Smoke test for Approach A: sample_reference + weight_at.
// Checks: (1) weight_at(pl_ref) recovers generate weight;
//         (2) mean weight_at(+1) ≈ generate(+1).weight (closure at fixed t).

#include "meradgen_molpol.hpp"
#include "meradgen_api.hpp"
#include "meradgen_globals.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <random>

using namespace meradgen;

static int g_fail = 0;

static void check(bool ok, const char* msg) {
  if (!ok) {
    std::cerr << "FAIL: " << msg << "\n";
    ++g_fail;
  } else {
    std::cout << "ok: " << msg << "\n";
  }
}

int main() {
  const double elab = 11.0;
  const double thetacm = 90.0 * pi / 180.0;
  const double phi = 0.0;
  const double pl_ref = 0.0;

  itest = 0;
  merad_init(elab);

  double vp[4];
  vpgen_from_angles(elab, thetacm, phi, vp);

  std::mt19937 rng(42u);
  std::uniform_real_distribution<double> uni(0.0, 1.0);

  // --- Self-consistency: weight_at(pl_ref) == sample weight (soft + hard) ---
  {
    bool saw_soft = false, saw_hard = false;
    for (int trial = 0; trial < 400 && !(saw_soft && saw_hard); ++trial) {
      double rand4[4] = {uni(rng), uni(rng), uni(rng), uni(rng)};
      MolPolEvent kin;
      if (!sample_reference(vp, rand4, kin, pl_ref))
        continue;

      WeightPieces w;
      if (!weight_at(pl_ref, vp, kin, pl_ref, w)) {
        check(false, "weight_at(pl_ref) ok");
        break;
      }

      const double rel = std::fabs(w.weight - kin.weight)
                         / std::max(std::fabs(kin.weight), 1e-30);
      if (kin.ich == 0 && !saw_soft) {
        std::cout << "  soft ich=0 kin.weight=" << kin.weight
                  << " w.weight=" << w.weight << " lr=" << w.lr << " rel=" << rel
                  << "\n";
        check(std::fabs(w.lr - 1.0) < 1e-9, "soft lr(pl_ref) ≈ 1");
        check(rel < 1e-9, "soft weight_at(pl_ref) matches sample");
        saw_soft = true;
      }
      if (kin.ich == 1 && !saw_hard) {
        std::cout << "  hard ich=1 kin.weight=" << kin.weight
                  << " w.weight=" << w.weight << " lr=" << w.lr << " rel=" << rel
                  << "\n";
        check(std::fabs(w.lr - 1.0) < 1e-9, "hard lr(pl_ref) ≈ 1");
        check(rel < 1e-9, "hard weight_at(pl_ref) matches sample");
        saw_hard = true;
      }
    }
    check(saw_soft, "saw soft event for self-check");
    check(saw_hard, "saw hard event for self-check");
  }

  // --- Closure: E_0[W(+1)] ≈ sitot(+1)/xs0(+1) ---
  {
    const int n = 200;
    int n_ok = 0;
    double sum_w = 0.0;
    double target = 0.0;
    bool have_target = false;

    for (int i = 0; i < n; ++i) {
      double rand4[4] = {uni(rng), uni(rng), uni(rng), uni(rng)};
      MolPolEvent kin;
      if (!sample_reference(vp, rand4, kin, pl_ref))
        continue;

      WeightPieces wp;
      if (!weight_at(+1.0, vp, kin, pl_ref, wp))
        continue;

      sum_w += wp.weight;
      ++n_ok;

      if (!have_target) {
        // One generate(+1) at same t gives the constant MERADGEN weight.
        MolPolEvent ev_p;
        double r2[4] = {uni(rng), uni(rng), uni(rng), uni(rng)};
        if (generate(+1.0, vp, r2, ev_p)) {
          target = ev_p.weight;
          have_target = true;
        }
      }
    }

    check(n_ok > n / 2, "enough accepted events for closure");
    check(have_target, "have generate(+1) target weight");

    const double mean = sum_w / static_cast<double>(n_ok);
    const double rel = std::fabs(mean - target) / std::max(std::fabs(target), 1e-30);
    std::cout << "  n_ok=" << n_ok << " mean W(+1)=" << mean
              << " target sitot/xs0=" << target << " rel=" << rel << "\n";
    // Loose: O(1/√N) statistical; 200 events → expect ~few %.
    check(rel < 0.15, "closure E_0[W(+1)] ≈ generate(+1).weight (15%)");
  }

  if (g_fail) {
    std::cerr << g_fail << " check(s) failed\n";
    return 1;
  }
  std::cout << "all checks passed\n";
  return 0;
}
