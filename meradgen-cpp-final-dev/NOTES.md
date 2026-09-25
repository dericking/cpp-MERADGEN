# meradgen-cpp-final-dev — NOTES

Working fork of production `meradgen-cpp-final/` for **Approach A**
(`sample_reference` / `weight_at`). Not the production library until
deliberately promoted.

Formerly named `meradgen-cpp-final-2`; renamed to `meradgen-cpp-final-dev`
for the `MERADGEN-CPP-MOLPOL-DEV` branch.

## Purpose

MolPol needs one shared radiative electron-pair sample and helicity
dependence only in weights. Production `meradgen()` both samples and
weights at the same `pl`, so `generate_pair(+1,−1)` can yield different
\((v,t_1,z)\) and different lab electrons.

This tree adds:

- `sample_reference(vp, rand4, kin, pl_ref=0)` — freeze kinematics / 4-vectors
- `weight_at(pl, vp, kin, pl_ref, WeightPieces&)` — soft/hard density LR
  without re-drawing CDFs or `vectrec`
- `DESIGN_WEIGHT_AT.md` — LR convention and caveats
- `tools/weight_at_smoke.cpp` — self-check + mean-weight closure

**Kernel physics sources** (`fsir.cpp`, `meradgen_main.cpp`, etc.) are
unchanged vs `meradgen-cpp-final`. Only the MolPol-facing wrapper
(`meradgen_molpol.*`) and CMake tests were extended.

## User / project decisions affecting this tree

1. Sample at \(P_{\mathrm{ref}}=0\) by default; reweight to \(\pm1\).
2. Recommended `WeightPieces::weight = LR * sitot_ref / xs0(P)` matches
   \(\mathbb{E}[W]=\mathrm{sitot}/\mathrm{xs0}\) (MERADGEN’s event weight).
3. For MolPol **analyzing power** on a shared track, use **`lr`** (or
   equivalent \(\propto\sigma_\pm\)), not `weight` — confirmed by
   `kinematics_examination` closure (\(A_{\mathrm{lr}}\approx-7/9\) at 90°).
4. Negative hard `fsir` left signed; `dens_ref≈0` → `weight_at` returns false
   (caller rerolls).
5. Closure / high-stats tests owned by `kinematics_examination/`, not here.
6. Consumers: `MolPolSandbox` CMake points here; later MolPol may FetchContent
   a tagged release once promoted.

## Build

```bash
cmake -S meradgen-cpp-final-dev -B meradgen-cpp-final-dev/build \
  -DMERADGEN_BUILD_TESTS=ON
cmake --build meradgen-cpp-final-dev/build -j
./meradgen-cpp-final-dev/build/weight_at_smoke
```

See `DESIGN_WEIGHT_AT.md`, `README_WORKING.md`.
