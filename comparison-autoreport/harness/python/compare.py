#!/usr/bin/env python3
"""Compare FORTRAN and C++ parity-driver event files. Exit 1 if physics bar fails."""
from __future__ import annotations

import argparse
import math
import re
import struct
import sys
from pathlib import Path

EVENT_RE = re.compile(r"^EVENT\s+(\d+)\s+START")
BRACE_RE = re.compile(r"\{([^}]*)\}")


def _parse_kv_block(path: Path, start: str, end: str) -> dict:
    """Parse a START/END key=value block. Missing block → empty dict."""
    header: dict = {}
    in_block = False
    if not path.is_file():
        return header
    for line in path.read_text().splitlines():
        s = line.strip()
        if s == start:
            in_block = True
            continue
        if s == end:
            break
        if not in_block or "=" not in s:
            continue
        key, val = s.split("=", 1)
        key, val = key.strip(), val.strip()
        if key in ("mode", "dump_policy"):
            header[key] = val
        elif key == "vpgen":
            header[key] = [float(x) for x in val.split()]
        elif key in ("max_rad", "max_calls", "calls", "radiative"):
            header[key] = int(float(val))
        else:
            header[key] = float(val)
    return header


def parse_header(path: Path) -> dict:
    """Kinematics and dump policy written at the top of a parity dump."""
    return _parse_kv_block(path, "HEADER START", "HEADER END")


def parse_footer(path: Path) -> dict:
    """Call / radiative counts written at the end of a parity dump."""
    return _parse_kv_block(path, "FOOTER START", "FOOTER END")


def parse_scalar(tok: str) -> float:
    t = tok.strip()
    if t.lower().startswith("0x"):
        u = int(t, 16)
        return struct.unpack("=f", struct.pack("=I", u))[0]
    return float(t)


def parse_list(inner: str) -> list[float]:
    return [parse_scalar(x) for x in inner.split()]


def parse_events(path: Path) -> list[dict]:
    events = []
    cur = None
    for line in path.read_text().splitlines():
        m = EVENT_RE.match(line)
        if m:
            if cur is not None:
                events.append(cur)
            cur = {"index": int(m.group(1)), "random": None, "vprad": None,
                   "phirad": None, "kin": None}
            continue
        if cur is None:
            continue
        if "RANDOM{" in line:
            inner = BRACE_RE.search(line)
            cur["random"] = [float(x) for x in inner.group(1).split()]
        elif "VPRAD{" in line:
            inner = BRACE_RE.search(line)
            cur["vprad"] = parse_list(inner.group(1))
        elif "PHIRAD{" in line:
            inner = BRACE_RE.search(line)
            cur["phirad"] = parse_list(inner.group(1))
        elif "KIN{" in line:
            inner = BRACE_RE.search(line)
            toks = inner.group(1).split()
            cur["kin"] = {
                "ich": int(float(toks[0])),
                "vgen": float(toks[1]),
                "t1gen": float(toks[2]),
                "zgen": float(toks[3]),
            }
    if cur is not None:
        events.append(cur)
    return events


def rel_delta(a: float, b: float) -> float:
    return abs(a - b) / max(abs(a), abs(b), 1e-30)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("fortran")
    ap.add_argument("cpp")
    ap.add_argument("--fail-rel", type=float, default=1e-6)
    ap.add_argument("--record-rel", type=float, default=1e-10)
    ap.add_argument("--bits", action="store_true",
                    help="official bar: fail if any VPRAD/PHIRAD float32 "
                         "component is not bitwise equal")
    args = ap.parse_args()

    fe = parse_events(Path(args.fortran))
    ce = parse_events(Path(args.cpp))
    if len(fe) != len(ce):
        print(f"FAIL event count fortran={len(fe)} cpp={len(ce)}", file=sys.stderr)
        return 1

    n_phys = 0
    n_rec = 0
    n_bit = 0
    n_kin = 0
    first_fail = None
    first_bits = None
    worst = (0.0, 0, "", 0)

    names = ["VPRAD[0]", "VPRAD[1]", "VPRAD[2]", "VPRAD[3]",
             "PHIRAD[0]", "PHIRAD[1]", "PHIRAD[2]", "PHIRAD[3]"]

    for f, c in zip(fe, ce):
        if f["index"] != c["index"]:
            print(f"FAIL index mismatch {f['index']} vs {c['index']}", file=sys.stderr)
            return 1
        if f["vprad"] is None or c["vprad"] is None or f["phirad"] is None or c["phirad"] is None:
            print(f"FAIL missing vectors event={f['index']}", file=sys.stderr)
            return 1
        fv = f["vprad"] + f["phirad"]
        cv = c["vprad"] + c["phirad"]
        event_phys = False
        event_rec = False
        for i, (a, b) in enumerate(zip(fv, cv)):
            if math.isnan(a) or math.isnan(b) or math.isinf(a) or math.isinf(b):
                if a != b:
                    event_phys = True
                    rd = float("inf")
                else:
                    continue
            else:
                rd = rel_delta(a, b)
            if rd > worst[0]:
                worst = (rd, f["index"], names[i], i)
            if rd > args.fail_rel:
                event_phys = True
                if first_fail is None:
                    first_fail = (f["index"], names[i], a, b, rd)
            elif rd > args.record_rel:
                event_rec = True
            if a != b and first_bits is None:
                first_bits = (f["index"], names[i], a, b, rd)
        if any(x != y for x, y in zip(fv, cv)):
            n_bit += 1
        if event_phys:
            n_phys += 1
        if event_rec:
            n_rec += 1
        if f.get("kin") and c.get("kin"):
            fk, ck = f["kin"], c["kin"]
            if fk["ich"] != ck["ich"] or any(
                fk[k] != ck[k] for k in ("vgen", "t1gen", "zgen")
            ):
                n_kin += 1

    print(f"events={len(fe)}")
    print(f"fail_rel={args.fail_rel} events_above_fail={n_phys}")
    print(f"record_rel={args.record_rel} events_above_record={n_rec}")
    print(f"events_with_any_bitwise_mismatch={n_bit}")
    if fe and fe[0].get("kin") is not None:
        print(f"events_with_kin_mismatch={n_kin}")
    if worst[1] == 0:
        print("worst_rel=0 (all compared components equal)")
    else:
        print(f"worst_rel={worst[0]:.6e} event={worst[1]} component={worst[2]}")
    if args.bits and first_bits:
        idx, name, a, b, rd = first_bits
        print(
            f"FAIL bits event={idx} {name} fortran={a:.16e} cpp={b:.16e} rel={rd:.6e}",
            file=sys.stderr,
        )
        return 1
    if first_fail and not args.bits:
        idx, name, a, b, rd = first_fail
        print(
            f"FAIL first event={idx} {name} fortran={a:.16e} cpp={b:.16e} rel={rd:.6e}",
            file=sys.stderr,
        )
        return 1
    print("PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
