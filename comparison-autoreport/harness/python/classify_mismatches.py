#!/usr/bin/env python3
"""Classify FORTRAN vs C++ event dumps (Stage 0a). Stdlib only. Does not run meradgen."""
from __future__ import annotations

import argparse
import struct
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from compare import parse_events, rel_delta  # noqa: E402

NAMES = ["VPRAD[0]", "VPRAD[1]", "VPRAD[2]", "VPRAD[3]",
         "PHIRAD[0]", "PHIRAD[1]", "PHIRAD[2]", "PHIRAD[3]"]


def f32_bits(x: float) -> int:
    return struct.unpack("=I", struct.pack("=f", x))[0]


def f32_hex(x: float) -> str:
    return f"0x{f32_bits(x):08X}"


def ulp32(a: float, b: float) -> int:
    """IEEE-754 binary32 ULP distance (finite values)."""
    ia = struct.unpack("=i", struct.pack("=f", a))[0]
    ib = struct.unpack("=i", struct.pack("=f", b))[0]
    if ia < 0:
        ia = 0x80000000 - ia
    if ib < 0:
        ib = 0x80000000 - ib
    return abs(ia - ib)


def bucket(rd: float) -> str:
    if rd == 0.0:
        return "0"
    if rd <= 1e-10:
        return "(0, 1e-10]"
    if rd <= 1e-8:
        return "(1e-10, 1e-8]"
    if rd <= 1e-6:
        return "(1e-8, 1e-6]"
    return "> 1e-6"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("fortran")
    ap.add_argument("cpp")
    ap.add_argument("--fail-rel", type=float, default=1e-6)
    ap.add_argument("--record-rel", type=float, default=1e-10)
    ap.add_argument("--top", type=int, default=20, help="print this many worst events")
    ap.add_argument("--all-bits", action="store_true",
                    help="print every bitwise-mismatch event (not only those above fail-rel)")
    args = ap.parse_args()

    fe = parse_events(Path(args.fortran))
    ce = parse_events(Path(args.cpp))
    if len(fe) != len(ce):
        print(f"event count fortran={len(fe)} cpp={len(ce)}", file=sys.stderr)
        return 1

    hist = {"0": 0, "(0, 1e-10]": 0, "(1e-10, 1e-8]": 0, "(1e-8, 1e-6]": 0, "> 1e-6": 0}
    n_bit = 0
    n_single_comp = 0
    n_rad = 0
    n_rad_bit = 0
    n_born_bit = 0
    ulp_hist: dict[int, int] = {}
    comp_hist = {n: 0 for n in NAMES}
    bit_events: list[tuple[int, bool, int, list]] = []
    worst: list[tuple[float, int, str, float, float, float]] = []
    fails = []

    for f, c in zip(fe, ce):
        fv = f["vprad"] + f["phirad"]
        cv = c["vprad"] + c["phirad"]
        radiative = any(x != 0.0 for x in f["phirad"])
        if radiative:
            n_rad += 1
        comps = []
        bit_comps = []
        event_worst = 0.0
        event_name = ""
        bit = False
        event_max_ulp = 0
        for i, (a, b) in enumerate(zip(fv, cv)):
            if a != b:
                bit = True
                u = ulp32(a, b)
                event_max_ulp = max(event_max_ulp, u)
                bit_comps.append((NAMES[i], a, b, rel_delta(a, b), abs(a - b), u))
                comp_hist[NAMES[i]] += 1
            rd = rel_delta(a, b)
            if rd > event_worst:
                event_worst = rd
                event_name = NAMES[i]
            if rd > args.record_rel:
                comps.append((NAMES[i], a, b, rd, abs(a - b)))
        hist[bucket(event_worst)] += 1
        if bit:
            n_bit += 1
            ulp_hist[event_max_ulp] = ulp_hist.get(event_max_ulp, 0) + 1
            bit_events.append((f["index"], radiative, event_max_ulp, bit_comps))
            if radiative:
                n_rad_bit += 1
            else:
                n_born_bit += 1
            if len(bit_comps) == 1:
                n_single_comp += 1
        if event_worst > args.fail_rel:
            fails.append((f["index"], event_name, event_worst, comps, radiative))
        if event_name:
            fi = NAMES.index(event_name)
            worst.append((event_worst, f["index"], event_name, fv[fi], cv[fi], abs(fv[fi] - cv[fi])))
        else:
            worst.append((0.0, f["index"], "", 0.0, 0.0, 0.0))

    worst.sort(reverse=True)
    print(f"events={len(fe)} radiative={n_rad}")
    print(f"bitwise_mismatch={n_bit} radiative={n_rad_bit} born={n_born_bit} single_component={n_single_comp}")
    print("max_ulp histogram (events):")
    for u in sorted(ulp_hist):
        print(f"  {u} ULP  {ulp_hist[u]}")
    print("mismatching components (count of components, not events):")
    for n in NAMES:
        if comp_hist[n]:
            print(f"  {n:10s} {comp_hist[n]}")
    print("max_rel histogram:")
    for k in ("0", "(0, 1e-10]", "(1e-10, 1e-8]", "(1e-8, 1e-6]", "> 1e-6"):
        print(f"  {k:16s} {hist[k]}")
    print(f"events_above_fail({args.fail_rel:g})={len(fails)}")
    for idx, name, rd, comps, rad in fails:
        print(f"FAIL event={idx} radiative={int(rad)} worst={name} rel={rd:.6e}")
        for n, a, b, r, ad in comps:
            if r > args.fail_rel:
                print(f"  {n} F={a:.16e} C={b:.16e} rel={r:.6e} abs={ad:.6e}")
    print(f"top {args.top} by max_rel:")
    for rd, idx, name, a, b, ad in worst[: args.top]:
        if rd == 0.0:
            break
        print(f"  event={idx} {name} rel={rd:.6e} abs={ad:.6e} F={a:.16e} C={b:.16e}")
    if args.all_bits:
        print(f"all {len(bit_events)} bitwise-mismatch events:")
        for idx, rad, max_u, comps in bit_events:
            print(f"EVENT {idx} radiative={int(rad)} max_ulp={max_u}")
            for n, a, b, r, ad, u in comps:
                print(
                    f"  {n} ulp={u} rel={r:.6e} abs={ad:.6e} "
                    f"F={a:.16e} {f32_hex(a)} C={b:.16e} {f32_hex(b)}"
                )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
