# MERADGEN Approach A — MolPolSandbox integration notes

Local sandbox only. Do **not** edit `_tmp/halla_molpol_sim/`.

## Goal

1. Sample radiation at \(P=0\); freeze the electron-pair lab kinematics.
2. Evaluate \(W_+\), \(W_-\) (and \(W_0\)) at that point via `weight_at`.
3. Shoot **one** coincidence electron pair into Geant4 (photon not tracked).
4. Record weights + MERADGEN kinematics on the ROOT tree.

## Files touched

| Path | Change |
|---|---|
| `CMakeLists.txt` | `add_subdirectory` → `../meradgen-cpp-final-2`; link `meradgen_cpp` |
| `include/MolPolMeradgen.hh` / `src/MolPolMeradgen.cc` | Approach A: `sample_reference` + `weight_at(±1)` |
| PGA / Event / IO / Messenger | `/MolPol/gen meradgen`, ROOT `evMerad*` branches |
| `macros/runexample.mac` | `gen meradgen`, 2.2 GeV, 500 events |

## How PGA calls MERADGEN

```text
GeneratePrimaries (gentype == "moller")
  → SampleMollerVertex (θcm, φ, XY/Z, MS, external eloss, Levchuk, zLum)
  → switch(RCtype):
       0: GenerateMollerBorn
       1: GenerateMollerAlexander  (structure-function path)
       2: GenerateMollerMeradgen   (sample_reference + weight_at lr±)
```

`fRadCorrFlag` removed — use `/MolPol/RCtype` only.
`/MolPol/radCorrections` remains as a legacy bool → RCtype 0/1.
## ROOT branches

| Branch | Content |
|---|---|
| `evMeradUsed` | 1 if MERADGEN path |
| `evMeradWeightsReady` | 1 when both `weight_at` calls succeeded |
| `evMeradW0` / `Wplus` / `Wminus` | Density LRs \(\mathrm{lr}_\pm=\sigma_\pm/\sigma_0\) (and average). Use these for \(A=(W_+-W_-)/(W_++W_-)\). **Not** `WeightPieces::weight` (`sitot/xs0`). |
| `evMeradIch`, `V`, `T1`, `Z` | frozen radiation kinematics |
| `evMeradE1`, `E2` | lab electron energies (GeV) |
| `evMeradRand4[4]` | sampling uniforms |

Also fills legacy `evUnpolWght`, `evPolPlus/MinusWghtZ` (scaled).

## Library

Uses `meradgen-cpp-final-2` (`sample_reference` / `weight_at`). See that tree’s
`DESIGN_WEIGHT_AT.md`.

## Build / run

```bash
cd MolPolSandbox/build   # or mkdir -p build && cd build && cmake ..
cmake --build . -j
./MolPol ../macros/runexample.mac
```
