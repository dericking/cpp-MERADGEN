# MERADGEN port status

This is the living status of the FORTRAN-to-C++ port. Update it when the
situation changes. The Cursor rules in `.cursor/rules/` do not copy this
file on purpose: a stale copy would be worse than none.

Other files have other jobs:

| File | Job |
|---|---|
| `ParityWork.md` | Checklist for making C++ match FORTRAN |
| `DevelopmentHistory.md` | How we found things (append-only; do not rewrite) |
| `TODO-FINAL.md` | Production port (`meradgen-cpp-final/`) |

Parity at physics 1e-6 is signed off. Two hex 1-ULP events are accepted
evaluation residuals, not a second formula.

**Last reviewed:** 2026-09-23

---

## Which tree to use

| If you need | Use | Why |
|---|---|---|
| The FORTRAN reference | `meradgen-fortran/` | Byte-identical to `meradgen10.tar` (checked 2026-08-31) |
| C++ that matches FORTRAN | `meradgen-cpp/` | Parity tree: FORTRAN literals, float `atan2`, `sngl` casts |
| C++ for MolPol / Geant4 | `meradgen-cpp-final/` | Production library (`add_subdirectory`). PDG constants, double `atan2`, no `sngl`. Compare this to the parity tree, not to FORTRAN |

`aj27` (without the debug prints) and float `atan2` are in `meradgen-cpp/`.
The old snapshot tree (`validation_checks/`) and the early `quick-check/`
smoke were deleted on 2026-09-23. The snapshot-only `aj31` term splits
were never copied over, so they are gone with that tree.

---

## Constants

Do not mix these in one tree.

| Tree | Constants |
|---|---|
| `meradgen-cpp/` | FORTRAN literals (this pass) |
| `meradgen-cpp-final/` | PDG 2024/25 (`TODO-FINAL.md`) |

---

## Signed off

- 50k shared stream, seed 20260831: **pass** at relative 1e-6.
- Official hex bar still fails two events, each by one float32 ULP. Both
  are evaluation, not physics (probe Stage 7). We accepted them. Do not
  start `aj31` to chase them.
  - **30363** `VPRAD[1]`: `sd` multiply association under a ~5e8:1
    cancellation (`z**2` vs left-associative `z*z`).
  - **32158** `VPRAD[3]`: inherits a 12-ULP `vgen` from `xsadd`.

---

## Known limitations

None of these blocks the parity sign-off. Most are already decided.

| Item | Status | Notes |
|---|---|---|
| Leftover `std::pow(x, 3\|4\|5)` in the parity `fsir.cpp` | Leave it | gfortran uses multiply trees; `libm` `pow` is ~1 ULP per site. That is `aj27` (already rewritten here) and the harmless `aj31` (`u**5`). It is not the 30363 fail. The production tree already switched those sites to `ipow3`/`ipow4`/`ipow5`. `pow(x,2)` is fine at `-O2`. |
| No default `CMAKE_BUILD_TYPE` | Real footgun | A bare `cmake -S . -B build` is `-O0`, and even `pow(x,2)` becomes a `libm` call. That build does not match the documented Release physics. One-line CMake default if we get tired of remembering `-DCMAKE_BUILD_TYPE=Release`. |
| `simps` when `b < a`, or `reps == aeps == 0` | Dead path | FORTRAN would integrate backwards / take a constant step. We do not. The only production call is `simpsx(1e-22, vmin, …)` with `vmin ≫ 1e-22`, so it never happens. Leave it (`ParityWork.md` §5). |
| `merad_init` is `double`; FORTRAN is `real` | Contract | Do not narrow inside the library. Parity drivers pass a float-rounded `elab`. `45` is already exact in float32; `10.6` and `10.7` are not — today’s 10.7 GeV run used the rounded HEADER value and was fine. |
| Not thread-safe | True of FORTRAN too | Namespace globals and `static` caches (`vmin`, `t1min`, `t1max`, `zmin`, the grids). Fine if MolPol generates one event at a time. Only matters if someone runs `meradgen` from several threads. |
| No CMake `install()` / `export()` | Packaging | The production tree already builds as a static library via `add_subdirectory`. `sirad` / `sinonr` are `sirad_out` / `sinonr_out` there. A proper install/export is for a later Geant4 packaging pass, not for matching FORTRAN. |

### Accepted residuals

None of these is a logic error.

| Event | What differs | Size | Why we kept it |
|---|---|---|---|
| 138994 | `PHIRAD[1]` | ~1.8e-12 abs (~0.06 ULP of float32) | Rounding tie after a 1–2 ULP libm difference |
| 30363 | `VPRAD[1]` | 1 float32 ULP | Same algebra, different association in `sd` (`z**2` vs `z*z`) under huge cancellation |
| 32158 | `VPRAD[3]` | 1 float32 ULP | Same class: `pow` vs `**`, Simpson `xsadd`, then `sngl` |

---

## Upstream oddities we reproduce on purpose

The port copies FORTRAN, including the weird bits.

- `fsir.f:310` has `2d0*v*(v-v)`, which is identically zero. It was
  probably meant to be `(v-u)`. We still emit zero.
- Three electron masses that do not agree with each other:

  | Symbol | Value | Where |
  |---|---|---|
  | `m` | `0.511000e-3` | main code |
  | `m2` | `0.261112e-6` | not equal to `m*m` |
  | `am2(1)` | `0.26110e-6` | `vacpol` |

- `merad_init` uses a `DATA` statement to initialise `COMMON` variables,
  which is not standard FORTRAN.

---

## How we test

| What | Where |
|---|---|
| FORTRAN ↔ live `meradgen-cpp/` | `comparison-autoreport/harness/` + `run.py` |
| Memo for any seed / N / kinematics | `comparison-autoreport/` |
| Formal tests of `meradgen-cpp-final/` | `tests/` (not started). Compare that tree to the parity tree, not to FORTRAN |

Parity-tree vs final-tree on smoke + hex2: four-vectors agree at relative
1e-6 … 2e-5, and `ich` matches (`meradgen-cpp-final/README.md`).

### 2026-08-31

Smoke is green at hex (all ten four-vectors match bit for bit). The 50k
stream (seed 20260831) is green at 1e-6 after we matched FORTRAN’s
default-real literals in C++ `vacpol`. The official hex bar still fails
two events, each by one float32 ULP (30363 `VPRAD[1]`, 32158 `VPRAD[3]`).
Those are evaluation residuals, not a second formula; we accepted them
and moved production work to `meradgen-cpp-final/`.

### 2026-09-15

Stepped up to 100k radiative events on the same shared stream (seed
20260831, 45 GeV) and wrote that up as a LaTeX memo in `comparison/`.
FORTRAN and C++ still agree at 1e-6, and the overlay histograms match
bin for bin. Six events differ by one ULP on a small `PHIRAD` component
— the same kind of leftover we already knew about. That memo is frozen.

### 2026-09-23

Turned the memo into something we can regenerate for any seed, event
count, or kinematics (`comparison-autoreport/`). The dumps now carry a
HEADER and FOOTER, so the report can be rebuilt from the run directory
without retyping the setup. A second 45 GeV stream (seed 20260923) had
seven residuals scattered through the run, so the cluster we saw on the
15th was just that seed. We also ran 100k at 10.7 GeV: still clean at
1e-6, thirty-one residuals, still one- or two-ULP `sngl` ties. Lower
energy makes `m/E` larger, so more values sit on a rounding boundary;
it is not a new disagreement. Deleted `validation_checks/` and
`quick-check/` — they were leftover.

**Saved this as my final parity report in the MERADGEN-CPP repo root.**
