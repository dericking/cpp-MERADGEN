# MERADGEN C++ (production)

This is the **production** tree. The FORTRAN was never edited to look
like C++. Matching FORTRAN was done in `meradgen-cpp/` (the parity
snapshot). This tree then walks C++ *away* from that snapshot on
purpose: PDG constants, double `atan2`, no `sngl`, and so on.

Compare this tree to `meradgen-cpp/`, not to `meradgen-fortran/`.
Number changes here are intentional physics/evaluation updates, not a
port regression.

Parity vs FORTRAN lives in `meradgen-cpp/` and is signed off at relative
1e-6 (two leftover hex 1-ULP events are evaluation residuals; see
`ParityWork.md` §8).

## What is live here (PDG 2024/25)

| Symbol | This tree | Parity tree (`meradgen-cpp/`) |
|---|---|---|
| `pi` | `acos(-1.0)` | `atan(1)*4` |
| `alfa` | `0.00729735256` | `0.729735e-2` |
| `m` | `0.51099895e-3` GeV | `0.511000e-3` |
| `m2` | `m * m` | `0.261112e-6` (not `m*m`) |
| `barn` | `0.38937966e6` | `0.389379e6` |
| `vacpol` `am2` | `{m2, mu2, tau2}` | FORTRAN DATA |
| `vacpol` loop | `10.0/9.0`, `ccc = 4.091` | FORTRAN `REAL*4` then promote |
| `vectrec` `atan2` | double | float then promote |
| `vprad`/`phirad` | `double` end-to-end | `sngl` (float cast) |
| `vpgen` | `double (px,py,pz,E)` | `float (px,py,pz,E)` |
| Z-grid `sd` | `-az*(z*z)-…` | left-associative `((-az)*z)*z` |
| `pow(x,3\|4\|5)` | `ipow3/4/5` multiply trees | `std::pow` (except `aj27`) |
| `aj47` | `v*(v-u)` (typo fixed) | `v*(v-v)` (zero; matches FORTRAN) |
| `aj28`/`aj29` log | `/2/m2**2` (typo fixed) | `/2/m2` (matches FORTRAN) |

## Parity-tree vs this tree (2026-08-31)

Same quads, `elab=45`, `thetacm=90°`, `phi=10°`, `pl=-1`. `ich` never
differs. Four-vector relative shifts are **~1e-6 … 2e-5**, not the 0.8%
τ-loop raw change (that is diluted inside `vacpol` then `sig`).

| Fixture | Events | Worst rel | Notes |
|---|---|---|---|
| smoke 10 | 10/10 differ | 1.03e-05 `PHIRAD[0]` (event 9, radiative) | Non-radiative sit at 1.08e-06 on `VPRAD[1]` (`m`/`pcm`) |
| hex2 (30362, 30363, 32157, 32158) | 4/4 differ | 2.36e-05 `PHIRAD[1]` (32158) | All radiative |

No event above 1e-4. That is “we updated physics on purpose,” not a port
regression.

## MolPol / Geant4

This tree is a **static library**, not a program. Geant4 compiles it as part
of MolPol. There is no standalone MERADGEN executable to ship.

In the MolPol `CMakeLists.txt`:

```cmake
add_subdirectory(path/to/meradgen-cpp-final meradgen_cpp)
target_link_libraries(MolPol PRIVATE meradgen_cpp)
```

Call `merad_init(elab)` once per beam energy. Draw four uniforms from
Geant4 (`G4UniformRand()` or the engine you already use) and pass them
as `rand4`. If `generate` / `generate_pair` returns `false`, reroll
those four numbers — do not clamp inside the kernel.

```cpp
#include "meradgen_molpol.hpp"

meradgen::merad_init(elab);  // GeV
double vp[4];
meradgen::vpgen_from_angles(elab, thetacm, phi, vp);  // angles in radians
// or fill vp as (px, py, pz, E) GeV CM yourself

double rand4[4] = { G4UniformRand(), G4UniformRand(),
                    G4UniformRand(), G4UniformRand() };
meradgen::MolPolEvent ev;
if (!meradgen::generate(pl, vp, rand4, ev)) {
  // reject / reroll
}
```

`include/meradgen_molpol.hpp`. Four-vectors are `(px, py, pz, E)` GeV CM.
`to_Epxpypz` / `to_pxpypzE` convert to the old MolPol `(E, px, py, pz)`
packing. Dual helicity with the same `rand4` is `generate_pair`.

`sirad` and `sinonr` are `sirad_out` / `sinonr_out` on the globals
(and on `MolPolEvent`). Born is `xs0`.

The library is not thread-safe (mutable globals). One event at a time per
thread, or serialize calls.

## Optional tools (not for Geant4)

```bash
cmake -S meradgen-cpp-final -B meradgen-cpp-final/build -DMERADGEN_BUILD_TOOLS=ON
cmake --build meradgen-cpp-final/build -j
```

That builds `meradgen_run` and `meradgen_quads` for local checks only.
