#!/usr/bin/env python3
"""Second-pass analyzer for closure_approach_a dumps.

Raw C++ dumps only; this script builds summaries / histograms.
See METHOD.md and meradgen-cpp-final-dev/DESIGN_WEIGHT_AT.md.
"""

from __future__ import annotations

import argparse
import math
import sys
from typing import Dict, List, Optional, Tuple


def parse_kv_line(line: str) -> Dict[str, str]:
    toks = line.split()
    out: Dict[str, str] = {}
    i = 0
    while i + 1 < len(toks):
        out[toks[i]] = toks[i + 1]
        i += 2
    return out


_META_KEYS = (
    "mode",
    "elab_GeV",
    "thetacm_deg",
    "phi_deg",
    "pl_sample",
    "pl_ref",
    "target_events",
    "seed",
    "radiative_only",
    "has_weight_at",
)


def load_dump(path: str) -> Tuple[Dict[str, str], List[Dict[str, str]]]:
    meta: Dict[str, str] = {}
    events: List[Dict[str, str]] = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                parts = line[1:].split()
                for i, tok in enumerate(parts):
                    if tok in _META_KEYS and i + 1 < len(parts):
                        meta[tok] = parts[i + 1]
                continue
            if line.startswith("event"):
                events.append(parse_kv_line(line))
    return meta, events


def fget(ev: Dict[str, str], key: str) -> Optional[float]:
    if key not in ev:
        return None
    try:
        return float(ev[key])
    except ValueError:
        return None


def mean_rms(xs: List[float], ws: Optional[List[float]] = None) -> Tuple[float, float]:
    if not xs:
        return float("nan"), float("nan")
    if ws is None:
        m = sum(xs) / len(xs)
        var = sum((x - m) ** 2 for x in xs) / len(xs)
        return m, math.sqrt(var)
    wsum = sum(ws)
    if wsum == 0.0:
        return float("nan"), float("nan")
    m = sum(x * w for x, w in zip(xs, ws)) / wsum
    var = sum(w * (x - m) ** 2 for x, w in zip(xs, ws)) / wsum
    return m, math.sqrt(max(0.0, var))


def histogram(
    xs: List[float],
    weights: Optional[List[float]],
    nbins: int,
    lo: float,
    hi: float,
) -> Tuple[List[float], List[float]]:
    if hi <= lo:
        hi = lo + 1.0
    counts = [0.0] * nbins
    for i, x in enumerate(xs):
        w = 1.0 if weights is None else weights[i]
        if x < lo or x >= hi:
            continue
        b = int((x - lo) / (hi - lo) * nbins)
        b = max(0, min(nbins - 1, b))
        counts[b] += w
    # normalize to probability mass
    s = sum(counts)
    if s > 0:
        counts = [c / s for c in counts]
    edges = [lo + (hi - lo) * i / nbins for i in range(nbins + 1)]
    return edges, counts


def l1_hist(a: List[float], b: List[float]) -> float:
    return sum(abs(x - y) for x, y in zip(a, b))


def analyze_paired(events: List[Dict[str, str]]) -> None:
    d_scat = [x for x in (fget(e, "dE_scat") for e in events) if x is not None]
    d_rec = [x for x in (fget(e, "dE_rec") for e in events) if x is not None]
    ep = [x for x in (fget(e, "Pp_E_scat") for e in events) if x is not None]
    e0 = [x for x in (fget(e, "P0_E_scat") for e in events) if x is not None]

    m_ds, r_ds = mean_rms(d_scat)
    m_dr, r_dr = mean_rms(d_rec)
    m_abs, _ = mean_rms([abs(x) for x in d_scat])

    print("=== paired bias diagnostic (same rand4; NOT Approach A) ===")
    print(f"N = {len(d_scat)}")
    print(f"ΔE_scat (Pp−P0) GeV: mean={m_ds:.6g}  rms={r_ds:.6g}  mean|Δ|={m_abs:.6g}")
    print(f"ΔE_rec  (Pp−P0) GeV: mean={m_dr:.6g}  rms={r_dr:.6g}")
    mp, rp = mean_rms(ep)
    m0, r0 = mean_rms(e0)
    print(f"E_scat Pp: mean={mp:.6g} rms={rp:.6g}")
    print(f"E_scat P0: mean={m0:.6g} rms={r0:.6g}")
    print("Nonzero ΔE ⇒ helicity-dependent sampling; Approach A freezes P=0 kin.")


def analyze_direct(events: List[Dict[str, str]], tag: str = "S") -> Dict[str, List[float]]:
    e_scat = [x for x in (fget(e, f"{tag}_E_scat") for e in events) if x is not None]
    e_rec = [x for x in (fget(e, f"{tag}_E_rec") for e in events) if x is not None]
    ms, rs = mean_rms(e_scat)
    mr, rr = mean_rms(e_rec)
    print(f"=== sample tag={tag} ===")
    print(f"N = {len(e_scat)}")
    print(f"E_scat GeV: mean={ms:.6g} rms={rs:.6g}")
    print(f"E_rec  GeV: mean={mr:.6g} rms={rr:.6g}")
    return {"E_scat": e_scat, "E_rec": e_rec}


def analyze_ref0(
    events: List[Dict[str, str]],
    has_wa: bool,
    direct: Optional[Dict[str, List[float]]] = None,
) -> None:
    analyze_direct(events, tag="P0")
    if not has_wa:
        print("=== closure reweight ===")
        print("SKIPPED: has_weight_at=0")
        return

    lrp: List[float] = []
    lrm: List[float] = []
    wp: List[float] = []
    wm: List[float] = []
    a_lr: List[float] = []
    a_xs0: List[float] = []
    e_scat: List[float] = []
    e_rec: List[float] = []
    n_wa_fail = 0

    for e in events:
        if fget(e, "wa_ok_pp") != 1.0 or fget(e, "wa_ok_pm") != 1.0:
            n_wa_fail += 1
            continue
        rp = fget(e, "Wp_lr")
        rm = fget(e, "Wm_lr")
        wplus = fget(e, "Wp_w")
        wminus = fget(e, "Wm_w")
        xp = fget(e, "Wp_xs0")
        xm = fget(e, "Wm_xs0")
        es = fget(e, "P0_E_scat")
        er = fget(e, "P0_E_rec")
        if None in (rp, rm, wplus, wminus, xp, xm, es, er):
            n_wa_fail += 1
            continue
        lrp.append(rp)
        lrm.append(rm)
        wp.append(wplus)
        wm.append(wminus)
        if rp + rm != 0.0:
            a_lr.append((rp - rm) / (rp + rm))
        if xp + xm != 0.0:
            a_xs0.append((xp - xm) / (xp + xm))
        e_scat.append(es)
        e_rec.append(er)

    print("=== Approach A reweight (P_ref=0 → ±1) ===")
    print(f"N with weight_at = {len(wp)}  (skipped/fail {n_wa_fail})")
    if not wp:
        return

    mwp, rwp = mean_rms(wp)
    mlrp, _ = mean_rms(lrp)
    print(f"mean Wp_w = {mwp:.6g}  rms={rwp:.6g}")
    print("  (all-channel E[Wp_w] should ≈ generate(+1).weight ≈ sitot/xs0;")
    print("   radiative-only is conditional on ich==1 and need not match)")
    print(f"mean Wp_lr = {mlrp:.6g}")

    if a_lr:
        ma, ra = mean_rms(a_lr)
        print(f"mean A_lr=(lr+−lr−)/(lr++lr−) = {ma:.6g}  rms={ra:.6g}")
        print("  (physical shared-track asymmetry from density LR)")
    if a_xs0:
        ma, ra = mean_rms(a_xs0)
        print(f"mean A_xs0=(xs0+−xs0−)/(xs0++xs0−) = {ma:.6g}  rms={ra:.6g}")
        print("  Born soft reference at 90° CM: −7/9 ≈ −0.777…")

    # Shape: E_ref[f LR] / E[LR]  vs direct P=+1
    ms_u, _ = mean_rms(e_scat)
    ms_w, _ = mean_rms(e_scat, lrp)
    mr_u, _ = mean_rms(e_rec)
    mr_w, _ = mean_rms(e_rec, lrp)
    print(f"E_scat: unweighted mean={ms_u:.6g}  LR-weighted={ms_w:.6g}")
    print(f"E_rec:  unweighted mean={mr_u:.6g}  LR-weighted={mr_w:.6g}")

    if direct and direct.get("E_scat"):
        d_scat = direct["E_scat"]
        d_rec = direct["E_rec"]
        md, _ = mean_rms(d_scat)
        mdr, _ = mean_rms(d_rec)
        print("=== vs direct P=+1 dump ===")
        print(f"direct E_scat mean={md:.6g}  Δ(LR-w − direct)={ms_w - md:.6g}")
        print(f"direct E_rec  mean={mdr:.6g}  Δ(LR-w − direct)={mr_w - mdr:.6g}")

        lo = min(min(e_scat), min(d_scat))
        hi = max(max(e_scat), max(d_scat)) + 1e-12
        _, h_rw = histogram(e_scat, lrp, 15, lo, hi)
        _, h_d = histogram(d_scat, None, 15, lo, hi)
        print(f"E_scat hist L1(reweight, direct) = {l1_hist(h_rw, h_d):.4f}  (0=identical)")
        print("Pass guide: |Δmean| ≲ few×σ/√N and L1 shrinking with N.")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("dump", help="closure_approach_a output file")
    ap.add_argument(
        "--direct-compare",
        default=None,
        help="direct-mode dump for spectrum closure vs LR-reweighted ref0",
    )
    args = ap.parse_args()

    meta, events = load_dump(args.dump)
    mode = meta.get("mode", "?")
    has_wa = meta.get("has_weight_at", "0") == "1"
    print(f"file={args.dump}")
    print(
        f"mode={mode} elab={meta.get('elab_GeV', '?')} "
        f"thetacm={meta.get('thetacm_deg', '?')} "
        f"has_weight_at={int(has_wa)} N={len(events)}"
    )

    if mode == "paired":
        analyze_paired(events)
    elif mode == "direct":
        analyze_direct(events, tag="S")
    elif mode == "ref0":
        direct = None
        if args.direct_compare:
            _, de = load_dump(args.direct_compare)
            direct = analyze_direct(de, tag="S")
        analyze_ref0(events, has_wa, direct)
    else:
        print(f"unknown mode in dump meta: {mode}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
