# stream_hex2 fixtures

Copied from `20260831_validation_stream/_scratch/quads.txt`
(seed 20260831, N=50000). The two events that fail the official hex bar
after the `vacpol` REAL*4 edit, plus their −1 neighbors.

## `stream_hex2_fails.txt` (isolated; driver events 1–2)

| Driver event | Stream event | Failed component |
|---|---|---|
| 1 | 30363 | VPRAD[1] 1 ULP |
| 2 | 32158 | VPRAD[3] 1 ULP |

## `stream_hex2_quads.txt` (isolated; driver events 1–4)

| Driver event | Stream event | Role |
|---|---|---|
| 1 | 30362 | neighbor (four-vectors match; Z-grid noise present) |
| 2 | 30363 | fail VPRAD[1] |
| 3 | 32157 | neighbor (KIN differs; four-vectors match) |
| 4 | 32158 | fail VPRAD[3] |

Regenerate without the fixture:

```
python3 validation_checks_new/harness/python/generate_stream.py \
  --seed 20260831 --n 50000 --only 30362,30363,32157,32158 -o ...
```
