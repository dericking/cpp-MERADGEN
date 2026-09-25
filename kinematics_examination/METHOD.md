# Approach A closure test

**Claim:** Sample radiative Møller kinematics once at \(P_{\mathrm{ref}}=0\),
freeze lab 4-vectors, and reweight to \(P=\pm1\) with `weight_at`. Lab
spectra (and mean asymmetry) must match sampling directly at \(P=+1\)
(and \(P=-1\)).

Harness: `closure_approach_a` + `analyze_closure.py`.
Library contract: `meradgen-cpp-final-dev/DESIGN_WEIGHT_AT.md`.

## Why

MERADGEN CDFs depend on `pl`. Same `rand4` at \(P=+1\) vs \(P=0\) draws
**different** \((v,t_1,z)\) — lab energies differ by
\(\mathcal{O}(100\,\mathrm{MeV})\) (`compare_helicity`,
`out/pol_compare_11GeV_rad5000.txt`).

Approach A shares **kinematics**, not randoms:

1. `sample_reference(vp, rand4, kin, pl_ref=0)`
2. Freeze lab momenta from `kin`
3. `weight_at(+1)` and `weight_at(-1)` at frozen \((v,t_1,z,\mathrm{ich})\)

## Library API (present in `meradgen-cpp-final-dev`)

```cpp
bool sample_reference(vp, rand4, kin, /*pl_ref=*/0);
bool weight_at(pl, vp, kin, pl_ref, WeightPieces& out);
```

`WeightPieces::lr` = density likelihood ratio \(\sigma_P/\sigma_{\mathrm{ref}}\).
`WeightPieces::weight` = \(\mathrm{LR}\cdot\mathrm{sitot}_{\mathrm{ref}}/\mathrm{xs0}(P)\)
with \(\mathbb{E}_{\mathrm{ref}}[W]=\mathrm{sitot}(P)/\mathrm{xs0}(P)\).

CMake sets `MERADGEN_HAS_WEIGHT_AT=1` when `weight_at` appears in the
headers. If linking an older tree without the API, the harness still
builds: `ref0` dumps kinematics only; `paired` / `direct` always work.

## Closure procedure

Equal-\(N\) dumps at fixed \(t\) (\(\theta_\mathrm{CM}=90^\circ\), \(\phi=0\)):

```bash
# Direct P=+1
./build/closure_approach_a --mode direct --pl 1 --radiative --n 1000 \
  --out out/closure_direct_plus.txt

# Reference P=0 + weight_at
./build/closure_approach_a --mode ref0 --radiative --n 1000 --seed 2 \
  --out out/closure_ref0.txt

python3 analyze_closure.py out/closure_ref0.txt \
  --direct-compare out/closure_direct_plus.txt
```

**Pass (spectra):** LR-weighted \(E_\mathrm{scat}\) / \(E_\mathrm{recoil}\)
from `ref0` agree with unweighted `direct` means/shapes within MC noise
(\(\mathbb{E}_{\mathrm{ref}}[f\cdot\mathrm{LR}]/\mathbb{E}[\mathrm{LR}]
=\mathbb{E}_{P=+1}[f]\)).

**Pass (weights):** mean `Wp_w` ≈ `generate(+1).weight` (constant at fixed
\(t\)) — same check as `weight_at_smoke`.

**Pass (asymmetry):** physical shared-track asymmetry uses the density LR,
\(A_{\mathrm{lr}}=({\mathrm{lr}}_+-{\mathrm{lr}}_-)/({\mathrm{lr}}_++{\mathrm{lr}}_-)\),
near soft Born \(A_{\mathrm{xs0}}\approx -7/9\) at \(90^\circ\) (hard events pull
slightly). Do **not** form \(A\) from `WeightPieces::weight` — that quantity
carries an extra \(1/\mathrm{xs0}(P)\) so \(\mathbb{E}[W]=\mathrm{sitot}/\mathrm{xs0}\).

## Modes

| Mode | Role | Needs `weight_at`? |
|---|---|---|
| `direct` | Sample at `--pl` | No |
| `ref0` | `sample_reference` + optional `weight_at(±1)` | For full closure |
| `paired` | Same `rand4` \(P=+\) vs \(P=0\) (negative control) | No |

Default kinematics: **11 GeV**; also `--elab 2.2` (MolPol optics).

## Paired bias (negative control)

`paired` proves same-`rand4` is **not** Approach A. Nonzero
\(\Delta E_\mathrm{lab}\) is why we freeze \(P=0\) kinematics.
