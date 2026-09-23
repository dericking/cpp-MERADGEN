# FORTRAN vs parity-Cpp auto-report

Generate the FORTRAN vs `meradgen-cpp/` (parity tree) memo for **any**
shared-stream simulation: seed, event count, and kinematics are arguments.

The parity drivers and stream tools live in `comparison-autoreport/harness/`
(CMake, inject, fixtures, `compare.py`).

The official compared objects are still the 32-bit `sngl` `vprad`/`phirad`
four-vectors. `meradgen-cpp-final/` is not part of this comparison.

C++ was changed to match FORTRAN. The FORTRAN reference is not edited
to look like the port.

Automated report cooked up by Cursor from original report put together. *All
report results should be reviewed carefully.* You can write your own scripts
to review the simulated data dumps but this simplifies the report writing part. 

_Note: It takes ~500K events to get ~100K radiative events. To run the FORTRAN 
and CPP simulations and write the dumps, and run the report script takes about 
~18 minutes on my workstation which is modest but by no means underpowered.

## One command

From the repo root:

```
python3 comparison-autoreport/run.py \
  --seed 20260831 \
  --n-radiative 100000 \
  --n-quads 500000 \
  --elab 45 --thetacm 90 --phi 10 --pl -1
```

That writes a self-contained run directory (default
`comparison-autoreport/_scratch/run/`) and compiles the PDF.

`--n-quads` is the number of meradgen **calls** in the stream. About 23% of
calls are radiative at the default kinematics, so the stream must be long
enough to yield `--n-radiative` events. Non-radiative calls still run
(SAVE'd state); they are not written as EVENT blocks.

## What the dump now carries

Each parity dump (`fortran_radiative.txt`, `cpp_radiative.txt`) starts with a
`HEADER` and ends with a `FOOTER`, so a report can be rebuilt from the files
without re-typing the run:

```
HEADER START
elab=...
thetacm=...
phi=...
pl=...
m=...
m2=...
ecm=...
pcm=...
vpgen=...
mode=full
max_rad=...
max_calls=...
dump_policy=radiative_only
HEADER END
EVENT ...
FOOTER START
calls=...
radiative=...
FOOTER END
```

`elab` is float32-rounded, matching the parity driver contract. Optional
driver arguments (used by `run.py`, ignored by old callers) are

```
driver_parity quads.txt output.txt [es14|full|hex] [max_rad] [max_calls] \
              [elab thetacm phi pl]
```

## Run directory

| File | Role |
|---|---|
| `config.json` | Command, seed, requested N, kinematics |
| `quads.txt` / `quads.meta` | Shared `rand4` stream (seed, n, sha256) |
| `fortran_radiative.txt` / `cpp_radiative.txt` | Radiative events; `EVENT` index is the 1-based stream call; HEADER/FOOTER as above |
| `toolchain.txt` | Compiler versions and CMake flags |
| `outliers.tsv` | Every event with a nonzero photon or electron delta |
| `summary.json` | All numbers the PDF interpolates (no hardcoded 100k / 427348 / event ids) |
| `overlay_radiative.png` / `delta_photon_electron.png` | Figures |
| `MERADGEN_<N>_radiative_FORTRAN_vs_CPP.pdf` | The memo |

## Report only, from an existing run

```
python3 comparison-autoreport/analyze.py --run-dir comparison-autoreport/_scratch/run
python3 comparison-autoreport/write_report.py --run-dir comparison-autoreport/_scratch/run
```

`analyze.py` prefers dump HEADER kinematics (`elab`, `m`, `m2`, `ecm`, `pcm`)
when reconstructing \(k_2\).

## Another kinematics example [MOLLER-energy central ray]

```
python3 comparison-autoreport/run.py \
  --seed 1 --n-radiative 200 --n-quads 2000 \
  --elab 10.6 --thetacm 90 --phi 0 --pl 1 \
  --out-dir comparison-autoreport/_scratch/elab10p6
```

`10.6` is not a float32-exact value; the driver stores the float-rounded
`elab` in the HEADER, and that is what the report prints.
