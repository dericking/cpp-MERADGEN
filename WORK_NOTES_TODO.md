# TODO — final C++ port

Punch list for the **production** C++ tree. Parity is signed off
(`ParityWork.md` §8): physics 1e-6 vs FORTRAN; two hex 1-ULP events
accepted as evaluation residuals.

Live tree: **`meradgen-cpp-final/`**. Keep **`meradgen-cpp/`** as the
parity snapshot. Compare final vs parity (not vs FORTRAN) so a number
change is “we updated physics on purpose,” not “the port went wacky.”

---

## Constants — revert FORTRAN literals to PDG

Files: `meradgen-cpp-final/src/cpp/globals.cpp`,
`meradgen-cpp-final/src/cpp/meradgen_core.cpp` (`vacpol`).

Uncomment PDG, comment FORTRAN. Do not mix the two sets.

| Symbol | Live now (FORTRAN, parity) | Restore (PDG 2024/25) |
|---|---|---|
| `pi` | `atan(1)*4` | `acos(-1.0)` is fine |
| `alfa` | `0.729735e-2` | `0.00729735256` |
| `m` | `0.511000e-3` GeV | `0.51099895e-3` GeV |
| `m2` | `0.261112e-6` (not `m*m`) | `m * m` |
| `barn` | `0.389379e6` | `0.38937966e6` |
| `vacpol` `am2[3]` | `0.26110e-6`, `0.111637e-1`, `3.18301` | `{m2, mu2, tau2}` |

`mu` / `tau` are already PDG and unused on the FORTRAN `vacpol` path. After
the `am2` flip they become live. Do not leave FORTRAN `am2` in place with PDG
`m`/`m2` — that is a mixed set.

- [x] `globals.cpp`: switch `pi`, `alfa`, `m`, `m2`, `barn`
- [x] `vacpol`: switch `am2` to `{m2, mu2, tau2}`
- [x] Re-run the parity-tree vs final-tree comparison (expect ~1e-6 relative
      from `m`/`m2`/`barn`/`alfa`, and ~0.8% on the τ loop in `vacpol`)
      **2026-08-31:** smoke 10 + hex2 4. `ich` matches. Four-vector rel
      1e-6 … 2e-5 (worst 2.36e-05 on 32158 `PHIRAD[1]`). τ-loop 0.8% is
      diluted; nothing above 1e-4 on `VPRAD`/`PHIRAD`. See the final README.
- [x] Record the new constants and the PDG year in the C++ README

---

## `vectrec` `atan2` — revert float to double

File: `meradgen-cpp-final/src/cpp/meradgen_main.cpp` (`vectrec`).

Uncomment double, comment float. FORTRAN does `phi = atan2(vpgen(2), vpgen(1))`
on `real` (float32). Double `atan2` after promoting the args is the more
accurate evaluation.

```cpp
// parity (live now):
// const float phi_f = std::atan2(vpgen[1], vpgen[0]);
// const double phi  = static_cast<double>(phi_f);

const double phi = std::atan2(static_cast<double>(vpgen[1]),
                              static_cast<double>(vpgen[0]));
```

- [x] Restore double `atan2` in `vectrec`
- [x] Compare final tree vs parity tree on a radiative event (expected
      ~1e-7 rad in `phi`, amplified in `VPRAD[0/1]` / `PHIRAD[1]`)
      Hex2 radiative events: worst 2.36e-05 on `PHIRAD[1]`. Combined with
      PDG/`sngl`/ipow, not an isolated atan2 measurement.

---

## Other parity choices to drop in the final tree

Not constants, but the same class of “match FORTRAN, then revert.”

- [x] Drop `sngl` / store `vprad`/`phirad` as `double` end-to-end
- [x] Remaining `std::pow(x, 3|4|5)` → multiply trees (accuracy, not just 138994)
      206 sites in `fsir.cpp` → `ipow3`/`ipow4`/`ipow5`. Fractional
      `**1.5`/`**2.5` and `pow(x,2)` stay (`pow(x,2)` folds at `-O2`).
- [x] `vacpol`: restore double literals (`10.0/9.0`, `ccc = 4.091`) in place of
      FORTRAN default `REAL*4` then promote (parity live in `meradgen-cpp/`)
- [x] CMake: static library for Geant4 (`add_subdirectory`). `-ffp-contract=off
      -fno-fast-math`. Does not force `CMAKE_BUILD_TYPE` when nested.
      Example binaries (`meradgen_run`, `meradgen_quads`) are
      `-DMERADGEN_BUILD_TOOLS=ON` only.
- [x] MolPol adapter (GeV, `(px,py,pz,E)`, injected `rand4`)
      `meradgen_molpol.hpp`: `generate` / `generate_pair`, `vpgen_from_angles`,
      `(E,px,py,pz)` packers. Kernel `meradgen()` now takes `double vpgen[4]`.
- [x] Z-grid `sd`: explicit `-az*(z*z)-…` (conventional grouping; parity tree
      left-associates `((-az)*z)*z`, which was the 30363 hex residual)

---

## Invalid kinematics — reject/reroll in MolPol, do not floor the kernel

The parity pass removed testing-era clamps in `vectrec` (`sqrt(max(al, 0))`,
`1e-100` floors on `denom`/`al1`) and the CDF `d > 1e-100` guards. Those
treated “numerically garbage” as “silently invent a 4-momentum.” They are
not in FORTRAN, they can fire on events FORTRAN still computed (C++ ULP
pushing `al8` slightly negative), and the denom floor could flip the sign
of `px`/`py`.

Do **not** put them back in a production `meradgen-cpp-*` tree.

The sound version of that idea: if `al1`/`al3`/`al8` are negative, a
denominator vanishes, a CDF bin width is ~0, or `vprad`/`phirad` is
NaN/Inf, **reject or reroll the event in the Geant4/MolPol layer**. The
generator kernel should keep producing the same (possibly NaN) result
FORTRAN would.

- [x] MolPol adapter: detect invalid reconstructed 4-momenta (NaN/Inf,
      negative radicands if exposed) and reject/reroll; do not clamp
      inside `vectrec` / CDF interpolation
      `event_is_finite` / `generate` return `false`. No kernel clamps.
- [x] Do not restore `sqrt(max(al, 0))`, `denom_safe`, or CDF `1e-100`
      guards in the kernel (verified: none in `meradgen-cpp-final/`)

---

## Optional, physics-review only

- [x] Upstream `fsir.f` `v*(v-v)` (identically zero; intended `v*(v-u)`).
      FORTRAN and the parity tree still emit zero. Production
      (`meradgen-cpp-final/`) uses `v*(v-u)` as of 2026-09-23.
- [x] Upstream `aj28`/`aj29` log `/2/m2` (intended `/2/m2**2`, same
      log as `aj21`/`aj22`). Production fixed 2026-09-23.
      See `README_MERADGEN_TYPO.md`.
