# Parity work list

Working notes for making `meradgen-cpp/` reproduce `meradgen-fortran/` event
output. This is a checklist, not a status board (`STATUS.md`) and not a
history log (`DevelopmentHistory.md`). Tick items as they are decided or done.

**Goal of this pass:** same events, same randoms, FORTRAN and C++ agree on
the `sngl` four-vectors at IEEE-754 **binary32** (hex). Relative 1e-6 is
still reported as a physics metric. `ES14.6` is a human dump.

**Later, not this pass:** a `meradgen-cpp-*/` tree that uses modern constants,
double `atan2`, explicit multiplies instead of `pow` only where we want
accuracy, no `sngl` casts, etc. That tree is compared to *this* parity tree,
not to FORTRAN, so we can tell “we changed physics on purpose” from “the
port went wacky.”

---

## Policy for this pass

1. Match FORTRAN behaviour, including its oddities (inconsistent electron
   masses, identically-zero `v*(v-v)` term, `REAL*4` `atan2` and `sngl`
   output casts).
2. Injected `rand4` stays. That is an interface change, not a physics change.
   The example driver (`run_main.cpp`) is not a parity driver.
3. Do not mix “make it match FORTRAN” and “make it better physics” in the
   same edit.
4. Dead comments / unused includes stay unless they get in the way.

---

## 1. Constants — use the FORTRAN literals

`globals.cpp` currently has PDG 2024/25 values. Every validated snapshot
used the FORTRAN set. For this pass, switch back.

| Symbol | FORTRAN (use these) | Current C++ (PDG) |
|---|---|---|
| `alfa` | `0.729735e-2` | `0.00729735256` |
| `m` | `0.511000e-3` | `0.51099895e-3` |
| `m2` | `0.261112e-6` (not `m*m`) | `m*m` |
| `barn` | `0.389379e6` | `0.38937966e6` |
| `vacpol` `am2[0,1,2]` | `0.26110e-6`, `0.111637e-1`, `3.18301` | `m2`, `mu2`, `tau2` |
| `pi` | `atan(1d0)*4d0` | was `acos(-1.0)` |

- [x] Restore FORTRAN literals in `globals.cpp` / `vacpol` (PDG lines commented beside them)
- [x] Keep PDG values commented so the later `meradgen-cpp-*` tree can pick them up in one place
- [x] `m2` stays the hardcoded FORTRAN number, not `m*m`
- [x] `pi` set to `atan(1)*4` to match `merad_init`

---

## 2. `vectrec` — float `atan2`

FORTRAN: `phi = atan2(vpgen(2), vpgen(1))` on `real` (float32).
Live C++ now matches that, then promotes to double for `sin`/`cos`. The
double `atan2` path is commented beside it for `TODO-FINAL.md`.

- [x] Restore the float `atan2`, then promote to double for `sin`/`cos`
- [x] Leave a comment that this is a parity choice; double path stays commented

---

## 3. `std::pow` vs FORTRAN `**`

Not “trivial precision” when it hits a `float` rounding boundary (event
138994 / `aj27`). The deep snapshot rewrote `aj27` with multiply trees;
that split is now live in canonical `fsir.cpp` (`ikey = 2`). Remaining
integer-3/4/5 `pow` sites (including `aj31`) stay as `pow` until a stream
comparison says we need them.

For this pass we do **not** have to rewrite all 206. Minimum that recovered
the known visible mismatch:

- [x] Back-port the `aj27` multiply/divide split from the (removed)
      spring 2026 snapshot **without** the `DBG` `std::cout` block. Live
      in `meradgen-cpp/src/cpp/fsir.cpp` (`ikey = 2`). Original `pow`
      line commented beside it.
- [x] `aj31` (`u**5`) and the remaining ~204 integer-3/4/5 sites: **deferred**
      until a stream comparison after this `aj27` backport. Campaign left
      `aj31` open (621/621 hex diffs, ~4e-14 %). If a visible event still
      diverges after `aj27`, do `aj31` here; otherwise later `meradgen-cpp-*`.

---

## 4. Defensive clamps (testing-era inserts)

None of these are in FORTRAN. They convert a NaN/Inf into a finite wrong
number.

| Site | C++ | FORTRAN |
|---|---|---|
| `vectrec` `sqrt(al1/al3/al8)` | `sqrt(max(al, 0))` | `sqrt(al)` |
| `vectrec` denom | `1e-100` floor | divide as written |
| V/T1/Z CDF interpolate | skip if `d > 1e-100` | divide unconditionally |
| `fspen` | `isnan`/`isinf` → `0` | no guard |

- [x] Remove `sqrt` / denom clamps in `vectrec`
- [x] Remove CDF `1e-100` guards (divide unconditionally, including
      negative `d`)
- [x] Remove `fspen` NaN/Inf fallback (and the dead `return 0.0`)
- [x] Keep the `fspen(x == 1.0)` *branch* fix (arithmetic-IF translation).
      `x <= 1.0` matches FORTRAN `if(x-1d0)4,4,5` (returns `f1`).
      Same for `x == 2.0`: `x <= 2.0` matches `if(x-2d0)6,6,7`.

---

## 5. Simpson integration

**Used.** FORTRAN `meradgen` calls

```fortran
call simpsx(1d-22, vmin, 10000, 1d-3, fsirv, xsadd)
```

to integrate additional soft bremsstrahlung into `sinonr`. `simpsx` is a
wrapper that calls `simps`. C++ has both (`meradgen_simpson.cpp`) and
`meradgen()` calls `simpsx` the same way. This is on the production path.

What is **not** called from `meradgen` itself: `simpt`, `simpu`, `simptx`,
`simpux`. Those are body-identical copies of `simps`/`simpsx` in
`meradgen10.f`. They are now in the C++ API as wrappers over `simps`/`simpsx`
so the port is complete; they may be removed later if we decide they are
clutter.

Known C++ `simps` translation issues (edge cases, not the xsadd call):

- [x] Port `simpt` / `simpu` / `simptx` / `simpux` (wrappers; FORTRAN bodies identical)
- [x] `b < a` and `reps == aeps == 0`: **do not arise** from
      `simpsx(1d-22, vmin, 10000, 1d-3, …)`. Documented and deferred.
- [x] Sentinel `1e16` vs FORTRAN `10.d16` (`1e17`): same call site; deferred.

**Why the edges never fire.** `simpsx` always calls
`simps(a=1e-22, b=vmin, reps=1e-3, aeps=1e-18)`.

- `vmin = 2 * Egmin * m = 0.02 * En * m`. At `En = 45` GeV that is
  ~`4.6e-4` ≫ `1e-22`. `b < a` would need `En ~ 10⁻¹⁷` GeV.
- FORTRAN `if (b-a) 1,2,1` integrates backwards when `b < a` and returns
  when `b == a`. C++ `if (b <= a) return` is therefore wrong in general
  and unused here (`b > a` always).
- `reps`/`aeps` are never both zero: `simpsx` hardcodes `1e-3` and `1e-18`.
  FORTRAN constant-step mode is unreachable from this wrapper.
- Sentinel: FORTRAN uses `10.d16` (`1e17`) with `==`; C++ uses `1e16` with
  `>=`. Only matters if an integrand sample equals the sentinel. Not worth
  hunting on `xsadd`. Revisit if a standalone `simps` caller appears.

No change to `meradgen_simpson.cpp` on this pass.

---

## 6. Other translation items (this pass)

- [x] `merad_init(elab)`: keep the C++ `double` signature. Do **not** narrow
      inside `merad_init`. Parity drivers pass a float-rounded energy
      (`static_cast<double>(45.f)`, or any `real` value FORTRAN would have
      stored). `45` is already float-exact; `10.6` is not.
- [x] `nn` vs `nn_t1` in the T1 loop: diagnostic counter only. Leave as-is.
- [x] `vprad`/`phirad` stored as `double` but filled with `float` casts.
      Matches FORTRAN `sngl` of the *values*. Fine for this pass.

---

## 7. How we will know we are done

Harness: `validation_checks_new/` (live trees, dated campaigns,
`VALIDATION_LOG.md`). Official bar: bitwise IEEE-754 binary32 of
`sngl` `VPRAD`/`PHIRAD` (`run.py` default `--prec hex`). Relative 1e-6
is still printed. 50k **calls** for the stream campaign.

- [x] `20260831_validation_smoke` green (10 smoke-fixture quads)
- [x] Stream campaign ran: seed 20260831, N=50000, fail at 1e-6.
      First run **FAIL** (3 events). After `vacpol` default-real literals
      in C++: **PASS** at 1e-6. 49999/50000 exact at `ES14.6`; leftover event
      32158 `VPRAD[3]` rel 5.71e-07.
- [x] Stream green at 1e-6. `vacpol` `10./9.` / `ccc=4.091` REAL*4
      promotion was the `xsvr` split. `L1f` matched. Remaining `pow` /
      `aj31` not required for this bar.
- [x] Official compare switched to hex (Stage 6 option A). Smoke **PASS**.
      Stream **FAIL** two 1-ULP events (30363 `VPRAD[1]`, 32158 `VPRAD[3]`).
      0 events above 1e-6. See `VALIDATION_LOG.md` and
      `20260831_validation_stream/fails.md`.
- [x] Those two hex fails localized: **evaluation**, not a second formula.
      30363 is `sd` multiply association (`z**2` vs `((-az)*z)*z`) under
      ~5e8:1 cancellation. 32158 inherits 12-ULP `vgen` from `xsadd`/`sitot`.
      Accepted. Do not rewrite remaining `pow` / `aj31` in this tree.

---

## 8. Stream residuals (2026-08-31)

Hex bar: 2 / 50000 events, each 1 float32 ULP on one small `VPRAD`
component. Same class as 138994 (tiny double noise, `sngl` tie). Physics
at rel 1e-6 is clean.

Localized 2026-08-31 (`stream_hex2_quads.txt` `--trace`):

| Event | What actually differs | Not |
|---|---|---|
| **30363** | Z-grid `sd`: FORTRAN `-az*z**2-…` vs C++ `((-az)*z)*z-…`. Inputs `az,bz,cz,z` match. 5.4e8:1 cancel → rel `1.29e-08` on every `aj*`/`sr*` as a common scale → `zgen` 973 ULPs → `sngl` 1 ULP on `VPRAD[1]` | Different `fsir` physics; `aj31` `u**5`; `vectrec`; `distsiv` |
| **32158** | `vgen` already 12 double ULPs from `xsadd` (6 ULP) + V-bin 58 `fsir.out` (1 ULP). T/Z inherit | A new kernel bug |

First chronological hex split on every event (including neighbors that
still match float32): `ikey=2` `aj31` 1 ULP (`pow` vs `**`). Harmless to
the V CDF. Neighbor 30362 has the same `sd` association with milder
cancel (~1e3:1) and still matches four-vectors.

**Parity pass is done.** Remaining integer `pow` stays in this tree.
Production work is `TODO-FINAL.md` / `meradgen-cpp-final/`, compared to
*this* tree, not to FORTRAN.

---

## Later tree
Parity achivement signed off 2026-08-31 (1e-6 green; two hex ULPs accepted as eval).
PDG constants, double `atan2`, drop `sngl`, restore double `vacpol`
literals: that list, in the final tree.

## Analysis on 10.7 GeV central ray radiative kinematics
100k radiative, seed 20260923, `elab=10.7` / `thetacm=90` / `phi=10` /
`pl=-1` (`comparison-autoreport/_scratch/elab10p7/`). Same stream as the
45 GeV seed-20260923 memo. `elab` is not float32-exact; the dump HEADER
stores `10.699999809265137`.

498308 calls, 20.1% radiative. **PASS** at rel 1e-6 (0 events above).
`ich` matches. Overlay bins identical. 31 events have a nonzero
photon or electron residual (28 `PHIRAD`, 3 `VPRAD`); all below the
physics bar. Same eval-residual class as 45 GeV, not a new kinematics
split.

### Explanation on increase in ULP differences 45GeV ==> 10.7GeV
The FORTRAN/C++ difference is still the usual few double ULPs from 
evaluation (pow vs **, multiply association, Simpson). Some of that 
noise comes from pieces that do not scale with beam energy. The electron 
mass is fixed, and `m/E` is about 4x larger at 10.7 GeV. Absolute eval 
noise therefore shrinks more slowly than the four-vectors. A double 
that sat inside a `float32` bin at 45 GeV sits closer to a `sngl` midpoint 
at 10.7 GeV, so more events flip by a single calculated bit