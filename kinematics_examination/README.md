# kinematics_examination

Inspection of the kinematics (and related weights / cross sections) that
MERADGEN produces for the \(P=+1\), \(P=0\), and \(P=-1\) states at fixed hard
process — and tools to validate **Approach A** (sample at \(P=0\), freeze
kinematics, reweight to \(P=\pm1\)).

User decisions and findings: [`NOTES.md`](NOTES.md).
Scientific contract: [`METHOD.md`](METHOD.md). Design context:
`MolPolMeradgenNEW/PLAN.md`.

| Tool | Role |
|---|---|
| `compare_helicity` | Same `rand4` at \(P=+1,0,-1\); dump \((v,t_1,z)\), weights, 4-vectors |
| `closure_approach_a` | Direct / ref0 / paired dumps with **lab** \(E\), \(\theta\) |
| `analyze_closure.py` | Histograms, \(\Delta E\) stats, reweight closure (when `weight_at` exists) |

Links **`meradgen-cpp-final-dev`** if present, else `meradgen-cpp-final`.
CMake sets `MERADGEN_HAS_WEIGHT_AT` when `weight_at` appears in the headers.

## Build

```bash
cmake -S kinematics_examination -B kinematics_examination/build
cmake --build kinematics_examination/build -j
```

Binaries:

- `kinematics_examination/build/compare_helicity`
- `kinematics_examination/build/closure_approach_a`

## Run — compare_helicity

```bash
./kinematics_examination/build/compare_helicity \
  --elab 11 --thetacm 90 --phi 0 \
  --n 200 --seed 1 --radiative \
  --out kinematics_examination/out/pol_compare_11GeV.txt
```

One dump: same quads, columns for \(P=+\), \(P=0\), \(P=-\).
Default angles are CM degrees. Quads from `--seed` unless `--quads FILE`.

Existing large dump: `out/pol_compare_11GeV_rad5000.txt`.

## Run — Approach A closure harness

Modes (see `METHOD.md`):

| `--mode` | Meaning |
|---|---|
| `direct` | Sample at `--pl` (default \(+1\)) |
| `ref0` | Sample at \(P=0\); call `weight_at(±1)` if compiled in |
| `paired` | Same `rand4` at \(P=+1\) and \(P=0\) (bias diagnostic) |

Default kinematics: **11 GeV**, \(\theta_\mathrm{CM}=90^\circ\), \(\phi=0\). Also use
`--elab 2.2` for MolPol optics energy.

```bash
# Smoke: paired bias diagnostic (works without weight_at)
./kinematics_examination/build/closure_approach_a \
  --mode paired --elab 11 --thetacm 90 --phi 0 \
  --n 500 --seed 1 --radiative \
  --out kinematics_examination/out/closure_paired_11GeV.txt

# Direct P=+1 sample
./kinematics_examination/build/closure_approach_a \
  --mode direct --pl 1 --elab 11 \
  --n 500 --seed 1 --radiative \
  --out kinematics_examination/out/closure_direct_plus_11GeV.txt

# P=0 reference (reweight columns appear only when weight_at is in the lib)
./kinematics_examination/build/closure_approach_a \
  --mode ref0 --elab 11 \
  --n 500 --seed 2 --radiative \
  --out kinematics_examination/out/closure_ref0_11GeV.txt

# 2.2 GeV optics energy
./kinematics_examination/build/closure_approach_a \
  --mode paired --elab 2.2 --n 200 --seed 1 --radiative \
  --out kinematics_examination/out/closure_paired_2p2GeV.txt
```

Analyze (Python 3, stdlib only):

```bash
python3 kinematics_examination/analyze_closure.py \
  kinematics_examination/out/closure_paired_11GeV.txt

# After weight_at lands:
python3 kinematics_examination/analyze_closure.py \
  kinematics_examination/out/closure_ref0_11GeV.txt \
  --direct-compare kinematics_examination/out/closure_direct_plus_11GeV.txt
```

## Status

- Links `meradgen-cpp-final-dev` when present. If that tree exports `weight_at`
  (see `DESIGN_WEIGHT_AT.md`), CMake sets `MERADGEN_HAS_WEIGHT_AT=1` and
  `ref0` dumps full reweight columns — run the analyzer with
  `--direct-compare` for spectrum closure.
- Without `weight_at`, `paired` / `direct` / kinematics-only `ref0` still run.
- Outputs under `out/` are gitignored.
