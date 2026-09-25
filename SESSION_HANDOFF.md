# Session handoff — MERADGEN ↔ MolPol (2026-09-24 → 2026-09-25)

**Give this file to the agent on the new workstation** (or open it and say
“continue from SESSION_HANDOFF.md”). Branch: `MERADGEN-CPP-MOLPOL-DEV` on
`https://github.com/dericking/cpp-MERADGEN.git`.

Also useful locally (this machine only): Cursor agent transcripts under
`~/.cursor/projects/.../agent-transcripts/` — not a substitute for pulling
this branch.

---

## 1. Goal of the session

Figure out how to use MERADGEN radiative corrections inside MolPol
(Geant4) so that:

- Coincidence **electron** kinematics are correct for spectrometer transport
- Helicity asymmetry is not biased by **different radiation** for \(P=+1\) vs \(P=-1\)
- The path is justifiable (closure / importance sampling)

Pristine MolPol clone `_tmp/halla_molpol_sim/` was **never modified**.

---

## 2. Core physics conclusions

### Same RNG ≠ same event across helicity

MERADGEN CDFs depend on `pl` (\(P = P_B P_T\)). Same `rand4` at \(P=+1\),
\(P=0\), \(P=-1\) produces **different** \((v,t_1,z)\), photons, and lab
electrons.

- 5000 radiative events @ 11 GeV, θ_CM=90°, φ=0: **0/5000** identical pairs
- Lab \(|\Delta E|\) scattered ~ **250 MeV** mean (Pp vs Pm) at 11 GeV; ~**40 MeV**
  at 2.2 GeV
- Angles barely move; **energies** matter for acceptance

### What is constant vs what varies (fixed `vpgen`)

| Quantity | Depends on \(P\)? | Varies event-by-event? |
|---|---|---|
| `xs0` = Born `sig(t,pl,0)` | yes | **no** (one value per \(P\)) |
| `sirad`, `sinonr`, MERADGEN `weight=sitot/xs0` | yes | **no** |
| `(v,t1,z)`, `phirad`, electron labs | yes (via CDF) | **yes** |

Born asymmetry at 90° CM: \(A_{xs0} \approx -0.7778\) vs ultra-rel limit
\(-7/9\). MERADGEN is exact finite-mass; \(-7/9\) is the UR limit.

### Approach A (chosen method)

1. Sample radiation once at **\(P_{\mathrm{ref}}=0\)**
2. Freeze lab electron pair (do **not** track the photon — none reach MolPol
   detectors)
3. Evaluate polarized densities at that frozen point → weights for \(P=\pm1\)
4. Shoot **one** coincidence pair into Geant4; ntuple carries \(W_\pm\)

**Do not** use `generate_pair` (Option B) for asymmetry — mismatched
electrons / acceptance.

### Weight convention (important)

Library `WeightPieces::weight = LR * sitot_ref / xs0(P)` matches
\(\mathbb{E}[W]=\mathrm{sitot}/\mathrm{xs0}\) (MERADGEN’s built-in weight).

For **shared-track analyzing power** use density **LR**:
\(A = (\mathrm{lr}_+ - \mathrm{lr}_-)/(\mathrm{lr}_+ + \mathrm{lr}_-)\).

Using `WeightPieces::weight` for \(A\) gives ~0.02 (RC-factor asymmetry),
**not** Born \(\approx-7/9\). Confirmed in closure and sandbox ROOT.

Absolute rate / luminosity weighting still to be reviewed carefully later
(user deferred).

---

## 3. Branch and trees

**Branch:** `MERADGEN-CPP-MOLPOL-DEV` (pushed to origin)

| Path | Role |
|---|---|
| `meradgen-cpp-final/` | Production library (unchanged kernel; still authoritative for MolPol releases until promoted) |
| `meradgen-cpp-final-dev/` | Working fork: `sample_reference` + `weight_at` (renamed from `*-final-2`) |
| `kinematics_examination/` | P+/P0/P− dumps + Approach A closure tools |
| `MolPolSandbox/` | Geant4 MolPol copy (no `.git`); RCtype integration |
| `_tmp/halla_molpol_sim/` | Pristine clone — **do not edit** |
| `MolPolMeradgenNEW/` | **Deleted** (superseded by sandbox adapter) |
| `MolPolMeradgenOLD/` | Gitignored prior glue |

Notes files (user decisions + status):

- `MolPolSandbox/NOTES.md`
- `MolPolSandbox/MERADGEN_INTEGRATION.md`
- `kinematics_examination/NOTES.md`
- `kinematics_examination/METHOD.md`
- `meradgen-cpp-final-dev/NOTES.md`
- `meradgen-cpp-final-dev/DESIGN_WEIGHT_AT.md`

---

## 4. Library API (`meradgen-cpp-final-dev`)

```cpp
bool sample_reference(vp, rand4, kin, /*pl_ref=*/0.0);
bool weight_at(pl, vp, kin, pl_ref, WeightPieces& out);
// WeightPieces: xs0, sinonr, fsir_z, dens, dens_ref, lr, sitot_ref, weight
```

- Soft LR: `sinonr(P)/sinonr(ref)`
- Hard LR: `fsir(...,P;ikey=0)/fsir(...,ref;ikey=0)`
- Kernel files (`fsir.cpp`, `meradgen_main.cpp`, …) **identical** to
  `meradgen-cpp-final`; only molpol wrapper + CMake tests added

Build smoke:

```bash
cmake -S meradgen-cpp-final-dev -B meradgen-cpp-final-dev/build -DMERADGEN_BUILD_TESTS=ON
cmake --build meradgen-cpp-final-dev/build -j
./meradgen-cpp-final-dev/build/weight_at_smoke
```

---

## 5. kinematics_examination

Purpose: **inspect** MERADGEN outputs at \(P=+1,0,-1\) and validate Approach A.

| Tool | Role |
|---|---|
| `compare_helicity` | Same `rand4`, dump all three \(P\) (raw; `--radiative` keeps ich=1) |
| `closure_approach_a` | `direct` / `ref0` / `paired` |
| `analyze_closure.py` | Histograms, \(A_{\mathrm{lr}}\), spectrum closure |

User dump prefs: raw only; φ=0; include P0 in same file; 5000 radiative @ 11 GeV
for stats. Outputs under `out/` gitignored.

Closure (N=500, all-channel): \(A_{\mathrm{lr}}\approx-0.77\); mean
`Wp_w` ≈ `generate(+1).weight` within ~0.4%.

```bash
cmake -S kinematics_examination -B kinematics_examination/build
cmake --build kinematics_examination/build -j
```

Links `meradgen-cpp-final-dev` when present.

---

## 6. MolPolSandbox

Copy of `halla_molpol_sim` without `.git`. Builds against
`../meradgen-cpp-final-dev`.

### UI

- `/MolPol/gen moller` (not a separate meradgen gentype)
- `/MolPol/RCtype 0|1|2` — **default 0**
  - 0 = Born
  - 1 = Alexander / Peking structure functions
  - 2 = MERADGEN Approach A
- Legacy `/MolPol/radCorrections` bool → RCtype 0/1 only (does not clear 2)

### PGA structure

```text
GeneratePrimaries (moller)
  → SampleMollerVertex (θcm, φ, vertex, MS, external eloss, Levchuk)
  → RCtype 0: GenerateMollerBorn
  → RCtype 1: GenerateMollerAlexander
  → RCtype 2: GenerateMollerMeradgen  (adapter: sample P=0, weight_at, LR→ntuple)
```

`fRadCorrFlag` **removed**. Born/Alexander share private
`GenerateMollerStructurePath(..., applyStructureRC)`.

### ROOT (`RCtype=2`)

`evMeradUsed`, `evMeradWeightsReady`, `evMeradW0/Wplus/Wminus` (**LR**),
`evMeradIch/V/T1/Z`, `evMeradE1/E2`, `evMeradRand4[4]`, plus legacy
`evUnpolWght` / `evPol±WghtZ` scaled from LR.

Photon not fired. Smoke @ 2.2 GeV: \(A_{\mathrm{lr}}\approx-0.772\).

```bash
cd MolPolSandbox/build   # or cmake -S MolPolSandbox -B MolPolSandbox/build
cmake --build . -j
./MolPol ../macros/runexample.mac   # RCtype 2, 2.2 GeV, 500 events
```

---

## 7. User decisions checklist (do not reverse casually)

1. Electrons matter; photons from MERADGEN are not detector objects.
2. Approach A (P0 sample + reweight), not dual sample.
3. RCtype integer; gen stays `moller`; no stacking 1+2.
4. Separate readable RC branches + shared vertex preamble.
5. Default RCtype = 0.
6. Justify method via closure (kinematics_examination); LR for \(A\).
7. Prototype in sandbox; pristine `_tmp/halla_molpol_sim` untouched; later
   PR to a **new MolPol branch**.
8. CMake pulls MERADGEN in-tree now; later FetchContent/GitHub **pin SHA/tag**.
9. Weight/luminosity convention to review later; LR for asymmetry for now.

---

## 8. Open / next work

1. **Weight semantics** — Document what MolPol should store for rates vs \(A\)
   (`lr` vs `lr*sitot_ref` vs `xs0*…`). Possibly add `xs0±` / `dens±` branches.
2. **Promote** `weight_at` from `meradgen-cpp-final-dev` → `meradgen-cpp-final`
   when ready.
3. **MolPol upstream PR** — Port sandbox PGA/IO/CMake to a new branch of
   `halla_molpol_sim` (not force onto pristine local clone).
4. **FetchContent** — Pin MERADGEN GitHub ref in MolPol CMake when public.
5. **Larger closure** — Raise N; optional thesis-quality spectrum L1 bounds.
6. **Remove per-event G4cout** RCtype prints (or gate behind verbose) before
   production runs.
7. Commit pending local tweak: `kinematics_examination/README.md` no longer
   points at deleted `MolPolMeradgenNEW/` (may be dirty vs remote).

---

## 9. New workstation bootstrap

```bash
git clone https://github.com/dericking/cpp-MERADGEN.git
cd cpp-MERADGEN
git checkout MERADGEN-CPP-MOLPOL-DEV
# read this file + MolPolSandbox/NOTES.md + kinematics_examination/NOTES.md
# Geant4 + ROOT needed for MolPolSandbox; library/kinematics tools are C++14 cmake
```

Prompt for the new agent:

> Continue MERADGEN–MolPol integration from `SESSION_HANDOFF.md` on branch
> `MERADGEN-CPP-MOLPOL-DEV`. Do not modify `_tmp/halla_molpol_sim` if present.
> Prefer `MolPolSandbox` + `meradgen-cpp-final-dev` + `kinematics_examination`.

---

## 10. Chat export note

There is no single “export this chat” button assumed here. Options:

1. **This file on the branch** (preferred handoff)
2. Local **agent transcripts** on the old machine (Cursor project
   `agent-transcripts/`) — copy if needed
3. Paste this summary into the new chat

End of handoff.
