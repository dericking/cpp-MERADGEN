# Development history

Laboratory notebook for the MERADGEN FORTRAN → C++ port. Entries are
append-only. Later findings that supersede earlier ones get a new dated
entry; old text is not rewritten.

Live status (which tree is authoritative, what is still open) lives in
`STATUS.md`, not here.

---

## 2026-03 through 2026-04 | first parity campaign (summary)

Work from roughly March–April 2026. This entry is a reconstruction written
on 2026-08-31 from the reports and snapshot trees under `validation_checks/`.
It is the starting record for this file, not a contemporaneous log.

### What the campaign was

MERADGEN is a FORTRAN 77 Møller-scattering radiative-correction event
generator (Afanasev, Ilyichev, Merenkov; CPC 2007 / hep-ph/0603027). The
port exists so MolPol (Geant4) can call it without a FORTRAN runtime. The
acceptance criterion used at the time was: given the *same* random numbers,
FORTRAN and C++ should produce the same events.

The FORTRAN owns its RNG (`URAND`, a linear congruential generator seeded
from `rnd.dat`). The C++ port was changed to take four injected uniforms
per event (`rand4[0..3]`) so both languages can consume a shared stream.
That injection is a validation interface, not a physics change.

Eight near-duplicate harness directories grew under `validation_checks/`
as the investigation deepened. They are frozen snapshots of that campaign,
not a test suite. The newest and most instrumented is
`validation_investigation_development_deep/` (directory mtime 2026-04-13).

### How events were compared

Two stream formats:

- Text quads (`r1 r2 r3 r4`, 10 decimal places) for development and
  single-event traces.
- Binary `MRADQF64` little-endian float64 streams for large-N runs
  (header + N×4 doubles). Typical campaign: consume until 50k *accepted
  radiative* events (`|PHIRAD| > 0`).

`r1` chooses the radiative vs non-radiative branch against `sitot`;
`r2` samples `t1`; `r3` samples `z`; `r4` chooses the `sl8` sign in
`vectrec`.

Important stream fact, learned the hard way: `meradgen` keeps `SAVE`d /
`static` grid state (`distsiv`, `distsit1`, `distsiz`, `sinonr`, `ikey`).
Event N in a stream is **not** reproducible by a cold call with only
quad N. Replay must run quads 1…N (`meradgen_replay_stream`). Event 138994
looked unreproducible until this was understood.

### Discovery 1 | float `atan2` in `vectrec`

Fortran:

```fortran
real vpgen(4)
phi=atan2(vpgen(2),vpgen(1))   ! REAL*4 atan2, then used as real*8
```

Early C++ promoted first:

```cpp
const double phi = std::atan2(static_cast<double>(vpgen[1]),
                              static_cast<double>(vpgen[0]));
```

The two `phi` values differ at ~1e-7 rad. That is then multiplied by large
kinematic coefficients (`sls*sl1*sl8`, `4*al3*al4 - s*al2*al7`) inside
`vectrec`, so a sub-microradian angle error becomes a ~1e-8 … 1e-7 absolute
error on `VPRAD[0/1]` and `PHIRAD[1]`. Seen on events 5 and 61.

**Fix used in the snapshots:** evaluate `atan2` in `float`, then widen.

```cpp
const float phi_f = std::atan2(vpgen[1], vpgen[0]);
const double phi  = static_cast<double>(phi_f);
```

This matches FORTRAN. It is a *parity* fix, not a physics improvement |
the float-width `atan2` is the less accurate of the two. The campaign
explicitly left a note that double `atan2` can be restored after parity
work is finished.

**Restart finding (2026-08-31):** `meradgen-cpp/` currently has the double
path active and the float path commented out. The validated snapshots still
have the float path. The fix was never (or was no longer) in the canonical
tree.

### Discovery 2 | `std::pow` vs Fortran `**` (aj27, event 138994)

Stream event 138994 was the first *visible* output mismatch after the
`atan2` fix: one `PHIRAD` component disagreed while a neighbouring event
(138990) did not. Tracing backwards, the first double-precision disagreement
was `aj27` (an `fsir` auxiliary, `ikey = 2` branch):

| | value |
|---|---|
| FORTRAN | `2.30801738251343450E+1` |
| C++ (`pow`) | `2.30801738251343380E+1` |

Last-bit of a `double`. Almost every other event absorbs this. 138994 does
not, because the subsequent `vectrec` quotient sits on a `float32` rounding
tie (`sngl` / `static_cast<float>` of `PHIRAD`). A 1 ULP shift upstream
flips which float32 bin the result lands in. The visible Δ is ~1.8e-12
absolute (~0.06 ULP of float32) on `PHIRAD[1]` only.

Root cause: the C++ translated every Fortran `**` as `std::pow`. Integer
exponents are *not* the same operation in the two compilers.

gfortran `-O2` expands them to multiply trees (read from assembly on
x86-64, GCC 11.4):

```
x**2  ->  x*x
x**3  ->  (x*x)*x
x**4  ->  t=x*x; t*t
x**5  ->  t2=x*x; t3=x*t2; t3*t2
```

g++ `-O2` only folds `pow(x,2)` to `x*x`. `pow(x,3)`, `pow(x,4)`, and
`pow(x,5)` remain calls to `libm` `pow`, which uses a different algorithm
(log/exp or a different multiply association) and can round 1 ULP away
from the multiply tree. At `-O0`, even `pow(x,2)` is a `libm` call.

**Worked example of compounding.** Suppose a kinematic factor `u` is used
as `u**5` in a denominator (this is exactly `aj31`) and also as `u**3` and
`u**4` in the numerator of neighbouring terms (`aj27` is built of those).
Each `pow` site can independently be 1 ULP off. Those terms are then added
into `sr1…sr10`, integrated by trapezoid over a 60- or 120-bin grid
(`distsiv` / `distsit1` / `distsiz`), and the CDF is inverted to sample
`vgen` / `t1gen` / `zgen`. A 1 ULP difference in a large `aj` value is
usually invisible in the integral. When the sampled point sits next to a
grid knot *and* a later `float` cast is within ~0.5 ULP of a bin edge, the
same 1 ULP becomes a different 4-momentum. That is 138994.

**Fix used in the deep snapshot:** split `aj27` into named intermediates
and replace `pow` with the gfortran multiply/divide pattern. After that,
441/441 `aj27` debug blocks matched FORTRAN at `.17E` for both 138994 and
the 138990 control.

**Restart finding (2026-08-31):** that fix exists only inside
`validation_investigation_development_deep/cpp_port/src/cpp/fsir.cpp`.
Canonical `meradgen-cpp/src/cpp/fsir.cpp` still has the monolithic `pow`
expression. A census of that file: 1,167 `std::pow` call sites; 206 of them
are integer exponents 3, 4 or 5 (the divergent class); 945 are `**2`
(benign at `-O2`); 16 are fractional (`**1.5`, `**2.5`), where gfortran
also calls `libm` so they should match.

### Discovery 3 | `aj31` / `u**5` (open at end of campaign)

Same method as `aj27`: split the monolithic `aj31` into intermediates
(`aj31_u_inner`, `aj31_log_coeff`, `aj31_sqrt_als`, `aj31_lin_comb`,
`aj31_log_num/den`, `aj31_dlog`, `aj31_log_term`, `aj31_u5 = u^5`).
Hex-compare FORTRAN vs C++ on event 138994.

Measured state in the April 13 reports:

- 621/621 occurrences: nonzero hexDiff on final `aj31`.
- 441/621 also differ on `aj31_u5` (`u**5`).
- Exact match on most log/linear pieces; persistent residual on
  `aj31` / `aj31_split` / `aj31_old` (200/441 nonzero, max ~4.1e-14 %).

Interpretation at the time: same class as `aj27` | the `u**5` powering
path | but not closed. There is no “441/441 match” conclusion analogous
to `aj27`. This was the last deep-session focus.

### Discovery 4 | leftover `PHIRAD[1]` float32 tie

Even after `aj27` matched, 138994 still differed on one returned float:

- FORTRAN `PHIRAD(2)` / C++ `phirad[1]`: −2.879630483e-05 vs −2.879630665e-05
- Abs Δ ≈ 1.8e-12

Diagnosed as 1–2 ULP in the double quotient in `vectrec` (libm `sin`/`cos`
and/or eval order) before the `sngl` cast, sitting on a float32 rounding
tie. `VPRAD` shares `sp`/`cp` and matched, so it was not a kinematic-input
mismatch. Treated as an accepted residual, not a logic bug. The
`validation_event_returns_updated_double/` tree promoted `vprad`/`phirad`
to `real*8`/`double` to study the same events without the float cast.

`vectrec` still casts outputs to `float` to match FORTRAN `real`. The
campaign left that as legacy behaviour, reason “unknown… left for now.”

### Discovery 5 | shared randoms, then large-N statistics

Before the shared stream, FORTRAN `URAND` and C++ `std::mt19937` were
independent generators, so event-by-event comparison was meaningless.
After sync:

- Quick-check (10 hardcoded quads, 3 radiative): relative agreement
  ~10⁻⁶ on `VPRAD`/`PHIRAD`.
- Large-N (50k accepted radiative): 99,999 / 100,000 events below
  10⁻¹⁰ percent; one event (the 138994 class) at 10⁻⁶ percent, traced
  as above.

Those figures were produced with *legacy* FORTRAN constants and with the
snapshot C++ (float `atan2`, and later the `aj27` multiply fix). They do
not describe canonical `meradgen-cpp/` as of 2026-08-31.

### Other notes from the campaign

- `itest` modes 1/2/3 short-circuit after sampling V / T1 / Z and were
  used as debug hooks; production is `itest = 0`.
- `ikey` / `SAVE`d distributions mean the first event after `merad_init`
  builds the V-grid; subsequent events with the same `t` reuse it. Dual
  helicity (`+pl` then `-pl` with the same `rand4`) therefore shares
  grids in a way a cold call does not.
- Output `weight` is `sngl(sitot/xs0)` | another float narrowing.
- Physical constants in the FORTRAN are internally inconsistent: `m`,
  `m2` (not equal to `m*m`), and `vacpol`’s `am2(1)` are three different
  electron masses. Snapshots reproduced all three. Canonical
  `meradgen-cpp/` later switched to PDG 2024/25 values (see next entry).

---

## 2026-06

Review of previous results and planning work for an actual parity
comparison to alleviate Jim's concerns about mismatches. 

- Start with clean FORTAN version
- Modify MERADGEN FORTRAN so that it uses the random-quad generator
- Generate a clean version of MERADGEN CPP which is as-close-as-possible replication of FORTRAN calculations
- 

---

## 2026-08-31 | restart: trees have drifted, FORTRAN is pristine

Restarted at end of summer prior to job switch. 

Compared `meradgen-fortran/` to a copy of the archived FORTRAN and 
to a fresh extraction of `meradgen10.tar`. All 17 files byte-identical 
(MD5 match). No unintended FORTRAN edits. The redundant archive directory 
was then removed.

Canonical `meradgen-cpp/` is **not** the tree that produced the campaign
numbers:

- Constants: PDG 2024/25 (`alfa`, `m`, `m2=m*m`, `barn`, lepton masses
  in `vacpol`) instead of the FORTRAN literals. Relative shifts ~1e-6 on
  `m`/`m2`/`barn`, ~3.5e-7 on `alfa`, ~0.8% on the τ loop in `vacpol`
  (`am2(3)` 3.18301 → τ² ≈ 3.157). This alone puts the canonical tree
  outside bitwise parity, and the τ change is a real physics change.
- `vectrec` `atan2`: double (campaign fix reverted).
- `aj27`: still the `pow` form (campaign fix never back-ported).
- `CMakeLists.txt` sets no default `CMAKE_BUILD_TYPE`, so a bare cmake
  configures `-O0` and even `pow(x,2)` stops matching gfortran.

Control-flow translation bugs found on this pass (not known to have been
fixed in any snapshot):

- `fspen(1.0)`: FORTRAN arithmetic `IF` sends `x==1` to the `f1` branch;
  C++ two-way tests send it into `log(1)*log(0)` → `NaN`.
- `simps` with `b < a`: FORTRAN integrates backwards; C++ returns 0.
- `simps` with `reps == aeps == 0`: FORTRAN constant-step mode; C++ has
  no such path.

`fsir.f:310` contains `2d0*v*(v-v)`, identically zero (likely a typo for
`(v-u)`). Both languages currently reproduce it.

Decision recorded: the port’s job is to reproduce the FORTRAN, including
its oddities, unless a later decision says otherwise. Constants (legacy
vs PDG) are undecided and block any new parity claim. See `STATUS.md`.

---

## 2026-08-31 | leftover hex fails are evaluation, not a second formula

Parity pass on live `meradgen-cpp/` vs `meradgen-fortran/` (FORTRAN
literals restored, float `atan2`, `aj27` multiply tree, `vacpol` default
`REAL*4` then promote). Official outputs are the eight `sngl` `VPRAD` /
`PHIRAD` components. 50k calls, seed 20260831, `elab=45.`.

Testing shows:

| Bar | Criterion | 50k result |
|---|---|---|
| rel 1e-6 | \|F−C\| / max(\|F\|,\|C\|) &lt; 10⁻⁶ | **PASS** (0 events above) |
| hex (official) | identical IEEE-754 binary32 bits | **FAIL** 2 events, 1 ULP each |

The two hex fails are 30363 `VPRAD[1]` (rel 7.5×10⁻⁸) and 32158
`VPRAD[3]` (rel 1.0×10⁻⁷). Both sit under 1e-6, so the physics bar never
sees them. Hex does, because 1 ULP of a small component is a different
bit pattern.

Those two events are FORTRAN vs C++ evaluation of the **same** algebra,
not a second cross section:

- **30363:** `vgen` / `t1gen` match. First Z-grid `fsir` call:
  `t, t1, v, z, az, bz, cz` bit-identical. Then every `aj*` / `sr*`
  differs by the same relative 1.291×10⁻⁸ (common scale = `sd`).
  FORTRAN `1d0/dsqrt(-az*z**2-…)` vs C++ `((-az)*z)*z`. The three
  ~10⁻⁶ terms cancel ~5.4×10⁸ : 1 (`disc ~ 8.2×10⁻¹⁵`). Python replay
  of the dumped coefficients reproduces each language’s `aj01` exactly.
  `zgen` moves 973 double ULPs → `sngl` lands in the next float32 bin.
  This is **not** the `ikey=2` `aj31 = π(…)/u**5` line (Z-grid `aj31`
  is `dz*dz*aj5`). Do not port the April `aj31` split.
- **32158:** `vgen` already 12 double ULPs from `xsadd` / `sitot`; T/Z
  inherit. Same class.

First chronological hex split on every traced event, including neighbors
that still match float32: `ikey=2` `aj31` at 1 double ULP (`pow(u,5)` vs
`u**5`). Harmless to the V CDF (`distsiv` identical).

Kinematics (`vgen` / `t1gen` / `zgen`) are dumped as `KIN` and are not
the pass/fail. ~6476 events differ there (worst ~10⁻¹⁰); 6474 of those
still match the float32 four-vectors.

Parity signed off at 1e-6. The two hex 1-ULP events are accepted
residuals. Remaining integer `pow` stays in `meradgen-cpp/`. Production
tree `meradgen-cpp-final/` started (PDG 2024/25, double `atan2`, no
`sngl`, double `vacpol` literals, explicit `-az*(z*z)` on `sd`). Compare
that tree to the parity snapshot, not to FORTRAN.
