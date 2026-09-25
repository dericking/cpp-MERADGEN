# MolPolSandbox — NOTES

Working Geant4/MolPol tree for MERADGEN Approach A integration.
**Not** the pristine clone: that lives at `_tmp/halla_molpol_sim/` (read-only;
do not modify). This sandbox was copied without `.git` so it is safe to edit
and to track on branch `MERADGEN-CPP-MOLPOL-DEV`. Build products under
`build/` are gitignored.

## Purpose

Prototype how MERADGEN (production C++ in `meradgen-cpp-final-dev/`) plugs
into MolPol’s Primary Generator and ROOT ntuple **without** changing
upstream `halla_molpol_sim` until a PR branch is ready.

## User decisions (recorded)

1. **Observable** — MolPol measures coincidence **electron pairs**, not
   photons. Do not track the MERADGEN photon as a Geant4 primary.
2. **Same radiation / same electrons** — Sample radiation once at \(P=0\),
   freeze the lab electron pair, evaluate helicity weights at \(P=\pm1\)
   (`weight_at`). Do **not** dual-sample with `generate_pair` for asymmetry.
3. **Generator type** — Stay on `/MolPol/gen moller`. RC is orthogonal.
4. **RCtype (integer)** — Prefer over a parallel bool:
   - `0` = none (Born kinematics; **default**)
   - `1` = Alexander / Peking structure-function approximation
   - `2` = MERADGEN Approach A
5. **Do not stack** RCtype 1 and 2.
6. **Layout** — Shared `SampleMollerVertex`, then separate readable branches
   `GenerateMollerBorn` / `GenerateMollerAlexander` / `GenerateMollerMeradgen`
   (like `gentype == "beam"`). Deprecate `fRadCorrFlag`; use `fRCtype` only.
   Legacy `/MolPol/radCorrections` maps to RCtype 0/1 only.
7. **CMake** — MolPol should pull MERADGEN into the build (`add_subdirectory`
   of `../meradgen-cpp-final-dev` for now). Later: FetchContent / GitHub
   pin by tag or SHA (not floating `main`).
8. **Weights for \(A\)** — Shared-track analyzing power uses density
   **likelihood ratios** (`lr±`), not MERADGEN’s `sitot/xs0`
   (`WeightPieces::weight`). Full weight convention to be reviewed later.
9. **Upstream push** — When ready, push sandbox-derived changes to a **new
   MolPol branch**, not directly onto pristine `_tmp/halla_molpol_sim`.

## What was implemented here

| Area | Status |
|---|---|
| CMake link to `meradgen-cpp-final-dev` | done |
| Adapter `MolPolMeradgen` (sample \(P=0\), `weight_at` ±1, LR → ntuple) | done |
| PGA `RCtype` 0/1/2 dispatch + three generators | done |
| ROOT `evMerad*` branches + Pol± Z weights from LR | done |
| Messenger `/MolPol/RCtype` | done |
| Smoke runs at 2.2 GeV (`macros/runexample.mac`, RCtype=2) | done |
| Photon gun | intentionally omitted |
| Proper absolute rate / luminosity weighting review | deferred |

See also `MERADGEN_INTEGRATION.md`, `SANDBOX.md`.

## Macro

`macros/runexample.mac`: `/MolPol/gen moller`, `/MolPol/RCtype 2`, 2.2 GeV
optics, short beamOn for sandbox tests.
