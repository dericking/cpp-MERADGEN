# MERADGEN C++ Port

This is the **parity** tree: C++ walked toward FORTRAN, not FORTRAN
walked toward C++. The reference is `meradgen-fortran/` and is never
edited. This tree reproduces it — FORTRAN literals, float `atan2`,
`sngl` casts, and the upstream oddities (identically-zero `v*(v-v)`,
inconsistent electron masses). Cleaner algebra here would be a physics
change. That work lives in `meradgen-cpp-final/`.

This directory is a standalone C++ build: a static library and a
runnable driver.

## Compiling and Running

### Compiling

From `meradgen-cpp/`:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

Artifacts produced:

- `build/libmeradgen_cpp.a`
- `build/meradgen_run`

### Running

From `meradgen-cpp/`:

```bash
./build/meradgen_run
```

The program prints event output to stdout as paired lines per event:

- `VPRAD` components
- `PHIRAD` components

## Changes Made

- Added injected random-quad flow to align event sampling with deterministic validation workflows:
  - `meradgen(...)` now accepts `rand4[4]`.
  - `vectrec(...)` now accepts `r4`.
  - Sampling decisions in V/T1/Z stages and final sign branch use injected random inputs.
- Updated `run_main.cpp` to generate random quads in-process using `std::mt19937` with a fixed seed (`12345`) and pass them into `meradgen(...)`.
  - This is for example only.
  - The parity harness (`comparison-autoreport/harness/`) generates a shared random-quad file that both FORTRAN and C++ consume.
- Kept the atan precision-parity change in `vectrec` (float-path `atan2` then promotion to double) in order to maintain consistency with legacy behavior used during validation.

Production (PDG constants, double `atan2`, no `sngl`) is `meradgen-cpp-final/`. Do not mix that work into this tree. 
