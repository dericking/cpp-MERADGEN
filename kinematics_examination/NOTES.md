# kinematics_examination — NOTES

Tools and dumps to **inspect MERADGEN kinematics and weights** for the
same hard process at polarization proxies \(P=+1\), \(P=0\), and \(P=-1\).
This directory is an examination / method-validation workbench, **not**
the Geant4 MolPol application.

## Purpose (user request)

Compare what MERADGEN produces for polarized vs unpolarized draws at fixed
lab energy and CM angles (default: 11 GeV, \(\theta_\mathrm{CM}=90^\circ\),
\(\phi=0^\circ\)), with an emphasis on:

- Whether the **same RNG quads** yield the **same radiation / electrons**
  across \(P\) (they do not).
- Cross sections / weights (`xs0`, `sirad`, `sinonr`, event `weight`) and
  radiation kinematics (`ich`, \(v\), \(t_1\), \(z\), `vprad`, `phirad`).
- Later: justify **Approach A** (sample at \(P=0\), freeze electrons,
  reweight with `weight_at` to \(\pm1\)).

Raw dumps only; calculations are a second pass (`analyze_closure.py`).

## User decisions recorded here

1. Dump **one file** for polarized (+/−) and include \(P=0\) in that same
   record so the \(P=0\) kinematics reference can be evaluated.
2. Prefer **raw data**, nothing precomputed in the dump (no match counts,
   no asymmetry columns in C++).
3. Fix \(\phi=0^\circ\) (scattering plane fixed).
4. Capture cross sections, weights, and kinematics (`ich`, \(v\), \(t_1\),
   \(z\), 4-vectors).
5. Run radiative (`ich==1` for all three \(P\)) samples: 50 then **5000**
   events for statistics.
6. Look at **emitted-photon** diffs/asymmetry (not `vpgen` — that is input).
7. Then look at **electron** lab kinematics (what MolPol transports);
   photons are not a detector concern.
8. Method: sample \(P=0\), freeze, reweight \(\pm\) — must be justified
   (closure / support / correct LR). Closure tools live here; library
   `weight_at` lives in `meradgen-cpp-final-dev/`.
9. For shared-track \(A\), use density **LR**, not `WeightPieces::weight`
   (`sitot/xs0`). Born \(A_{xs0}\approx-7/9\) at 90° CM (ultra-relativistic
   limit); MERADGEN `xs0` is exact finite-mass.

## Tools

| Binary / script | Role |
|---|---|
| `compare_helicity` | Same `rand4` at \(P=+1,0,-1\); one dump file |
| `closure_approach_a` | `direct` / `ref0` / `paired` modes for Approach A |
| `analyze_closure.py` | Second-pass histograms, \(A_{\mathrm{lr}}\), spectrum closure |

Build links **`meradgen-cpp-final-dev`** when present (else
`meradgen-cpp-final`). See `README.md`, `METHOD.md`.

## Key findings (summary)

- Same `rand4` ⇒ **0/5000** identical electron pairs or photons across \(P\).
- Lab \(\Delta E\) scattered \(\sim\) hundreds of MeV at 11 GeV (∼tens of MeV
  at 2.2 GeV) — enough to matter for spectrometer acceptance.
- At fixed `vpgen`, `xs0` / `sirad` / `sinonr` / MERADGEN `weight` are
  **constant per \(P\)** (Born / integrated totals), not per-event
  differentials.
- Born \(A_{xs0}\approx-0.7778\) vs \(-7/9\); RC-inclusive \(A_{\mathrm{sitot}}\)
  similar.
- Approach A closure: mean \(W(+1)\) matches `generate(+1).weight` at the
  percent level; \(A_{\mathrm{lr}}\approx-0.77\) near Born.

Outputs under `out/` are gitignored dumps.
