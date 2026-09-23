# stream_fails3 fixtures

Copied from `20260831_validation_stream/_scratch/quads.txt`
(seed 20260831, N=50000, `random.random()`, `{:.16e}`).

## `stream_fails3_quads.txt` (isolated; driver events 1–3)

| Driver event | Stream event | Failed component |
|---|---|---|
| 1 | 22992 | PHIRAD[2] |
| 2 | 34660 | VPRAD[1] |
| 3 | 40504 | PHIRAD[0] |

## `stream_fails3_neighbors.txt` (isolated; driver events 1–9)

| Driver event | Stream event | Role |
|---|---|---|
| 1 | 22991 | neighbor |
| 2 | 22992 | fail |
| 3 | 22993 | neighbor |
| 4 | 34659 | neighbor |
| 5 | 34660 | fail |
| 6 | 34661 | neighbor |
| 7 | 40503 | neighbor |
| 8 | 40504 | fail |
| 9 | 40505 | neighbor |

Regenerate without the fixture:

```
python3 validation_checks_new/harness/python/generate_stream.py \
  --seed 20260831 --n 50000 --only 22992,34660,40504 -o ...
```
