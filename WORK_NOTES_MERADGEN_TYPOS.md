# MERADGEN typos in `fsir.f`

Upstream oddities found in `meradgen-fortran/fsir.f`. The FORTRAN
reference is not edited. The parity tree (`meradgen-cpp/`) keeps the
upstream forms. The production tree (`meradgen-cpp-final/`) corrects
them. Compare production to parity, not to FORTRAN.

---

## 1. `aj47`: `v*(v-v)` → `v*(v-u)`

`ikey=2` formula for `aj47`, FORTRAN line 310. The coefficient of
`aj45` contains

```
2d0*v*(v-v)
```

That product is identically zero. Intended factor: `v*(v-u)`.

`(v-u)` is what you get after substituting `s-t-4m²`. In this generator
`u` is defined on line 25 as

```
u = v - s - t + 4 m²
```

so `v-u = s+t-4m²` (Mandelstam leftover once the photon takes `v`).
The same `aj47` block already uses a live `2d0*s*(v-u)` on line 311.

| Tree | `aj47` factor |
|---|---|
| `meradgen-fortran/` | `v*(v-v)` (zero) |
| `meradgen-cpp/` | `v*(v-v)` (zero) |
| `meradgen-cpp-final/` | `v*(v-u)` (fixed 2026-09-23) |

---

## 2. `aj28` / `aj29`: log denominator `/2/m2` → `/2/m2**2`

Same `ikey=2` block. The product log

```
dlog( (s-2m2)*(s-2m2+sqrt(als)) / 2 / m2?  - 1 )
```

appears four times:

| Site | FORTRAN line | Denominator | Verdict |
|---|---|---|---|
| `aj21` | 204 | `2 m2**2` | OK |
| `aj22` | 207 | `2 m2**2` | OK |
| `aj28` | 231 | `2 m2` | typo |
| `aj29` | 238 | `2 m2` | typo |

With \(s\) and \(m^2\) both GeV², the product is GeV⁴; only `/m2**2`
makes the argument of `log` dimensionless (and matches `aj21`/`aj22`).

Production C++ (`meradgen-cpp-final/src/cpp/fsir.cpp` lines 221–222):
`/2.0/std::pow(m2,2)`. Parity and FORTRAN keep `/2/m2`.
