# MERADGEN Project Workspace

Workspace containing the legacy FORTRAN MERADGEN source, a C++ port, and
FORTRAN vs C++ comparison reports.

**C++ was changed to match FORTRAN, not the other way around.**
`meradgen-fortran/` is the untouched reference (byte-identical to
`meradgen10.tar`). The parity tree (`meradgen-cpp/`) was walked toward
that FORTRAN — literals, float `atan2`, `sngl` casts, even upstream
oddities. The production tree (`meradgen-cpp-final/`) then walks C++
*away* from FORTRAN on purpose (PDG constants, double precision, and
selected typo fixes). We never edit the FORTRAN to look like the port.

## Top-level layout

| Path | Role |
|---|---|
| `meradgen10.tar` | Original upstream archive |
| `meradgen-fortran/` | FORTRAN reference (do not edit) |
| `meradgen-cpp/` | Parity C++ (match FORTRAN) |
| `meradgen-cpp-final/` | Production C++ for MolPol/Geant4 |
| `comparison-autoreport/` | Shared-stream memos (`run.py`) + parity `harness/` |
| `technote_*.pdf` | Citable notes (tracked; see `.gitignore`) |
| `WORK_NOTES_*.md` | Status, parity checklist, development log, typos, TODO |
| `MolPolMeradgenOLD/` | Old MolPol glue — local only (gitignored) |

## What we dumped

These are **gone** from the workspace (do not recreate the forked-tree pattern):

| Removed | What it was |
|---|---|
| `validation_checks/` | Spring 2026 snapshot forks (~eight full copies) |
| `quick-check/` | Early 10-event smoke |
| `validation_checks_new/` | Dated parity campaigns + scratch (2026-09-23) |
| `comparison/` | Frozen one-off 100k memo scripts (superseded) |

What replaced them:

- Parity drivers / inject / fixtures → `comparison-autoreport/harness/`
- Any-seed / any-kinematics memos → `comparison-autoreport/run.py`
- Scratch dumps → `comparison-autoreport/_scratch/` (gitignored)

## Notes

| File | Role |
|---|---|
| `WORK_NOTES_STATUS.md` | Live status and which tree is authoritative |
| `WORK_NOTES_PARITY.md` | Checklist for matching FORTRAN |
| `WORK_NOTES_DEVELOPMENT.md` | Append-only discovery log |
| `WORK_NOTES_TODO.md` | Production punch list |
| `WORK_NOTES_MERADGEN_TYPOS.md` | Upstream typos; fixed only in `meradgen-cpp-final/` |
| `technote_meradgen_port_compare_45GeV.pdf` | 100k radiative memo, 45 GeV |
| `technote_meradgen_port_comparison_10pt7GeV.pdf` | 100k radiative memo, 10.7 GeV |
| `technote_meradgen_fixes.pdf` | Production-side typo / fix note |

## Building the parity C++ port

```bash
cmake -S meradgen-cpp -B meradgen-cpp/build -DCMAKE_BUILD_TYPE=Release
cmake --build meradgen-cpp/build -j
```

Production library: see `meradgen-cpp-final/README.md` (`add_subdirectory` into MolPol).

## Running a comparison memo

```bash
python3 comparison-autoreport/run.py \
  --seed 20260831 \
  --n-radiative 100000 \
  --n-quads 500000 \
  --elab 45 --thetacm 90 --phi 10 --pl -1
```

That builds the harness, runs FORTRAN and parity C++ on one shared stream,
and writes a PDF under `comparison-autoreport/_scratch/`. Details:
`comparison-autoreport/README.md`.

Formal tests of `meradgen-cpp-final/` go under `tests/` later.
