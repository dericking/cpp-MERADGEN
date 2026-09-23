#!/usr/bin/env python3
"""Find the first hex-field mismatch in Stage 3 grid traces. Stdlib only."""
from __future__ import annotations

import argparse
import struct
import sys
from pathlib import Path


def hex_to_f64(h: str) -> float:
    t = h.strip()
    if t.lower().startswith("0x"):
        t = t[2:]
    u = int(t, 16)
    return struct.unpack("=d", struct.pack("=Q", u))[0]


def tokens(line: str) -> list[str]:
    return line.split()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("fortran")
    ap.add_argument("cpp")
    args = ap.parse_args()
    fl = Path(args.fortran).read_text().splitlines()
    cl = Path(args.cpp).read_text().splitlines()
    n = min(len(fl), len(cl))
    if len(fl) != len(cl):
        print(f"line_count fortran={len(fl)} cpp={len(cl)} (comparing first {n})")
    event = 0
    n_mis = 0
    first = None
    for i in range(n):
        a, b = fl[i], cl[i]
        if a.startswith("EVENT"):
            event = int(a.split()[1]) if a.split()[1:] else event
        if a == b:
            continue
        ta, tb = tokens(a), tokens(b)
        n_mis += 1
        rel = None
        if len(ta) >= 3 and ta[0] == "D" and len(tb) >= 3 and tb[0] == "D":
            try:
                fa, fb = hex_to_f64(ta[2]), hex_to_f64(tb[2])
                rel = abs(fa - fb) / max(abs(fa), abs(fb), 1e-30)
            except (ValueError, struct.error):
                rel = None
        rec = (i + 1, event, a, b, rel)
        if first is None:
            first = rec
            print(f"FIRST line={i+1} event={event}")
            print(f"  F {a}")
            print(f"  C {b}")
            if rel is not None:
                print(f"  rel={rel:.6e}")
    print(f"events_seen={event} differing_lines={n_mis}")
    if first is None:
        print("TRACE_MATCH")
        return 0
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
